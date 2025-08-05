import argparse
import re
import os
import shutil
import random

import numpy as np
import prismtoolbox as ptb
from pathlib import Path
from valis import valtils, preprocessing, registration
from utils_registration import customValis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform registration using Valis algorithm")
    parser.add_argument("--directory_IHC_reg", type=str, 
                        help="Path to the directory containing the IHC slides processed for registration")
    parser.add_argument("--directory_HE_reg", type=str, 
                    help="Path to the directory containing the H&E slides processed for registration")
    parser.add_argument("--contours_IHC_dir_reg", type=str,
                        help="Path to the directory containing the contours extracted from IHC slides and processed for registration")
    parser.add_argument("--contours_HE_dir_reg", type=str,
                        help="Path to the directory containing the contours extracted from H&E slides and processed for registration")
    parser.add_argument("--tmp_dir", type=str, default="./tmp",
                        help="Path to the directory where the temporary files will be saved")
    parser.add_argument("--result_dir", type=str,
                        help="Path to the directory where the results will be saved")
    parser.add_argument("--slide_ext", type=str, default="ome.tif",
                        help="Extension of the slide files processed for registration")
    parser.add_argument("--reference_slide", type=str, default="IHC", choices=["HE", "IHC"],
                        help="Reference slide for the registration")
    args = parser.parse_args()
    
    slides_IHC = os.listdir(args.directory_IHC_reg)
    slides_HE = os.listdir(args.directory_HE_reg)
    
    micro_reg_fraction = 0.25
    
    for slide_IHC in slides_IHC:
        tdx_nb = slide_IHC.split("-B0")[0]
        print(f"Processing Tdx Number {tdx_nb}")
        slide_HE = slides_HE[np.where([slide_IHC.split("-B0")[0] in s for s in slides_HE])[0].item()]
        
        # Load the slides and contours
        slide_IHC_path = os.path.join(args.directory_IHC_reg, slide_IHC)
        slide_HE_path = os.path.join(args.directory_HE_reg, slide_HE)
        contours_IHC_path = os.path.join(args.contours_IHC_dir_reg, f"{slide_IHC.split(f'.{args.slide_ext}')[0]}.pkl")
        contours_HE_path = os.path.join(args.contours_HE_dir_reg, f"{slide_HE.split(f'.{args.slide_ext}')[0]}.pkl")
        
        if os.path.exists(args.tmp_dir):
            shutil.rmtree(args.tmp_dir)
        Path(args.tmp_dir).mkdir(parents=True, exist_ok=True)
        
        # Create symbolic links to the slides and WSI objects with tissue contours
        WSI_object_dict = {}
        Path(os.path.join(args.tmp_dir, "slides")).mkdir(parents=True, exist_ok=True)
        for slide_type, slide_path, contours_path in zip(["HE", "IHC"], [slide_HE_path, slide_IHC_path],
                                                        [contours_HE_path, contours_IHC_path]):
            os.symlink(slide_path, os.path.join(args.tmp_dir, "slides", f"slide_{slide_type}.{args.slide_ext}"))
            WSI_object = ptb.WSI(slide_path, engine="tiffslide")
            WSI_object.load_tissue_contours(contours_path)
            WSI_object_dict[slide_type] = WSI_object
        
        # Create a custom Valis processor to create the mask using the tissue contours from the WSI object
        class IcarusProcessor(preprocessing.StainFlattener):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.WSI_object_dict = WSI_object_dict
                
            def create_mask(self):
                slide_type = re.search(rf"slide_(\w+).{args.slide_ext}", self.src_f).group(1)
                WSI_object = self.WSI_object_dict[slide_type]
                mask = 255 * np.array(WSI_object.visualize(self.level, black_white=True))
                return mask
        
        # Perform simple rigid registration to find the best seed
        errors = []
        for seed in range(1, 4):
            random.seed(seed)
            np.random.seed(seed)
                        
            with valtils.HiddenPrints():
                registrar = customValis(args.tmp_dir, "./tmp/reg",
                                        img_list=[f"./tmp/slides/slide_HE.{args.slide_ext}", f"./tmp/slides/slide_IHC.{args.slide_ext}"],
                                        crop="overlap",
                                        reference_img_f=f"./tmp/slides/slide_{args.reference_slide}.{args.slide_ext}",
                                        align_to_reference=True,
                                        create_masks=True,
                                        name=f"{slide_IHC.split('-B0')[0]}_{seed}",
                                        check_for_reflections=False,
                                        max_processed_image_dim_px=500,
                                        max_non_rigid_registration_dim_px=2000,
                                        non_rigid_registrar_cls=None)
                                            
                
                _, _, error_df = registrar.register(brightfield_processing_cls=IcarusProcessor,
                                                    brightfield_processing_kwargs={"adaptive_eq": True, "with_mask":True})
                errors.append(error_df["rigid_rTRE"].dropna().values[0])
            print(f"Error: {error_df['rigid_rTRE'].dropna().values[0]}")

        seed = int(np.argmin(errors) + 1)
        random.seed(seed)
        np.random.seed(seed)
        
        # Perform the registration
        with valtils.HiddenPrints():
            registrar = customValis(args.tmp_dir, args.result_dir,
                                    img_list=[f"./tmp/slides/slide_HE.{args.slide_ext}", f"./tmp/slides/slide_IHC.{args.slide_ext}"],
                                    crop="overlap",
                                    reference_img_f=f"./tmp/slides/slide_{args.reference_slide}.{args.slide_ext}",
                                    align_to_reference=True,
                                    create_masks=True,
                                    name=slide_IHC.split('-B0')[0],
                                    check_for_reflections=False,
                                    max_processed_image_dim_px=500,
                                    max_non_rigid_registration_dim_px=2000)
            
            rigid_registrar, non_rigid_registrar, error_df = registrar.register(brightfield_processing_cls=IcarusProcessor,
                                                                                brightfield_processing_kwargs={"adaptive_eq": True, 
                                                                                                               "with_mask":True})
            
            img_dims = np.array([slide_obj.slide_dimensions_wh[0] for slide_obj in registrar.slide_dict.values()])
            min_max_size = np.min([np.max(d) for d in img_dims])
            img_areas = [np.multiply(*d) for d in img_dims]
            max_img_w, max_img_h = tuple(img_dims[np.argmax(img_areas)])
            micro_reg_size = np.floor(min_max_size*micro_reg_fraction).astype(int)
            
            try:
                # Perform high resolution non-rigid registration using 25% full resolution
                micro_reg, micro_error = registrar.register_micro(max_non_rigid_registration_dim_px=micro_reg_size,
                                                                  reference_img_f=f"./tmp/slides/slide_{args.reference_slide}.{args.slide_ext}",
                                                                  align_to_reference=True)
            except Exception as e:
                print(f"Error: {e}, skipping micro registration")
                micro_reg, micro_error = None, None
        shutil.rmtree(args.tmp_dir)
    registration.kill_jvm()
        