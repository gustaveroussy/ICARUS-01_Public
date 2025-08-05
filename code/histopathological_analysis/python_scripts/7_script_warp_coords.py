import argparse
import os

import shapely
import numpy as np
import prismtoolbox as ptb
from valis import registration
from utils_registration import customValis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract patches from slides with contours already extracted")
    parser.add_argument("--directory_IHC_reg", type=str, help="Path to the directory containing the IHC slides processed for registration")
    parser.add_argument("--directory_HE_reg", type=str, help="Path to the directory containing the HE slides processed for registration")
    parser.add_argument("--contours_target_dir_reg", type=str,
                        help="Path to the directory containing the contours of the target slides processed for registration")
    parser.add_argument("--registrar_dir", type=str, help="Path to the directory containing the results of the registration")
    parser.add_argument("--patches_data_dir", type=str, help="Path to the directory containing the patches to warp")
    parser.add_argument("--results_dir", type=str, help="Path to the directory where the warped patches will be saved")
    parser.add_argument("--engine", type=str, default="tiffslide", help="Engine to use to read the slides")
    parser.add_argument("--slide_ext", type=str, default="ome.tif", help="Extension of the slides")
    parser.add_argument("--target_modality", type=str, default="IHC", choices=["HE", "IHC"],
                        help="Target slide for the registration")
    parser.add_argument("--file_format", type=str, default="h5", choices=["geojson", "h5"],
                        help="Format of the files to save the patches")
    parser.add_argument("--bbox_origin_dir", type=str, default=None, 
                        help="Path to the directory containing the bounding boxes of the tissue of interest used for the registration"
                        "(use when the patches were extracted prior to the registration)")
    args = parser.parse_args()
    origin_modality = "HE" if args.target_modality == "IHC" else "IHC"
    dir_origin = args.directory_HE_reg if origin_modality == "HE" else args.directory_IHC_reg
    dir_target = args.directory_HE_reg if args.target_modality == "HE" else args.directory_IHC_reg
    slides_target = os.listdir(dir_target)
    for slide_origin in os.listdir(dir_origin):
        print(f"Processing {slide_origin}")
        WSI_object_origin = ptb.WSI(os.path.join(dir_origin, slide_origin), engine=args.engine)
        
        if os.path.exists(os.path.join(args.results_dir, f"{WSI_object_origin.slide_name}.h5")):
            continue
        
        WSI_object_origin.load_patches(os.path.join(args.patches_data_dir, f"{WSI_object_origin.slide_name}.h5"))
        
        if args.bbox_origin_dir is not None:
            coords = shapely.MultiPoint(WSI_object_origin.coords)
            bbox = ptb.utils.data_utils.read_json_with_geopandas(os.path.join(args.bbox_origin_dir, f"{WSI_object_origin.slide_name}.geojson")).iloc[0].geometry
            intersect = coords & bbox
            if intersect.is_empty:
                raise ValueError(f"No patches in the bouding box for the slide {WSI_object_origin.slide_name}")
            else:
                WSI_object_origin.coords = np.stack([geom.coords[:][0] for geom in shapely.affinity.translate(intersect, xoff=-bbox.bounds[0], yoff=-bbox.bounds[1]).geoms]).astype("int64")
                
        tdx_nb = WSI_object_origin.slide_name.split("-B0")[0]
        
        registrar_f = os.path.join(args.registrar_dir, tdx_nb, "data", f"{tdx_nb}_registrar.pickle")
        registrar = registration.load_registrar(registrar_f)
        
        matching_target_slide = slides_target[np.where([slide_origin.split("-B0")[0] in s for s in slides_target])[0].item()]

        source_slide = registrar.get_slide(f"./tmp/slides/slide_{origin_modality}.{args.slide_ext}")
        target_slide = registrar.get_slide(f"./tmp/slides/slide_{args.target_modality}.{args.slide_ext}")
        warped_coords = source_slide.warp_xy_from_to(WSI_object_origin.coords, target_slide).round().astype("int64")
    
        WSI_object_target = ptb.WSI(os.path.join(dir_target, matching_target_slide), engine="tiffslide")
        
        WSI_object_target.load_tissue_contours(os.path.join(args.contours_target_dir_reg, f"{WSI_object_target.slide_name}.pkl"))
        valid_coords = []
        valid_indices = []
        patch_level = WSI_object_origin.coords_attrs["patch_level"]
        patch_size = WSI_object_origin.coords_attrs["patch_size"]

        valid_coords = []
        valid_indices = []
        unique_indices_set = set()
        for cont in WSI_object_target.tissue_contours:
            coords_in_cont, indices_in_cont =  WSI_object_target.extract_patches_roi(patch_level, patch_size,
                                                                                contour=cont,
                                                                                contours_mode="four_pt_hard",
                                                                                coord_candidates=warped_coords,
                                                                                rgb_threshs=(0, 255),
                                                                                percentages=(1, 1),
                                                                                return_indices=True)
            indices_in_cont = np.array(indices_in_cont)
            indices_to_keep = [k for k in range(len(indices_in_cont)) if indices_in_cont[k] not in unique_indices_set]
            valid_coords.extend(coords_in_cont[indices_to_keep])
            valid_indices.extend(indices_in_cont[indices_to_keep])
            for idx in indices_in_cont[indices_to_keep]:
                unique_indices_set.add(idx)
        
        
        assert len(valid_indices) == len(np.unique(valid_indices))
        WSI_object_target.coords = np.array(valid_coords)
        WSI_object_target.coords_attrs = WSI_object_origin.coords_attrs
        
        WSI_object_target.save_patches(args.results_dir, file_format=args.file_format)
        
        print(f"Number of patches warped that are valid: {len(valid_coords)}")
        print(f"Percentage of patches warped that are valid: {len(valid_coords)/len(WSI_object_origin.coords)}")
        
        WSI_object_origin.coords = WSI_object_origin.coords[valid_indices] # be careful not sorted
        WSI_object_origin.save_patches(args.results_dir, file_format=args.file_format)

registration.kill_jvm() # Kill the JVM