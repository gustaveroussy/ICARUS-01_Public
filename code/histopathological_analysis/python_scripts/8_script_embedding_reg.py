import argparse
import sys
import os
import logging

import pandas as pd
import numpy as np
import prismtoolbox.wsiemb as ptb_emb
from prismtoolbox.utils.data_utils import read_h5_file
from pathlib import Path

logging.basicConfig(stream=sys.stdout, level=logging.INFO)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract embeddings from WSI")
    parser.add_argument("--directory_HE_reg", type=str, help="Path to the directory containing the H&E slides processed for registration")
    parser.add_argument("--directory_IHC_reg", type=str, help="Path to the directory containing the IHC slides processed for registration")
    parser.add_argument("--directory_patches", type=str, 
                        help="Path to the directory containing the original patches and the warped patches (script 7)")
    parser.add_argument("--directory_cell_detection_reg", type=str, 
                        help="Path to the directory containing the cell detection processed for registration (script 5)")
    parser.add_argument("--results_dir", type=str, help="Path to the directory where the results will be saved")
    parser.add_argument("--tdx_to_remove", type=str, default=None, help="TDX to remove")
    parser.add_argument("--slide_ext", type=str, default="ome.tif", help="Extension of the slides")
    parser.add_argument("--patch_size", type=int, default=128, help="Size of the patches")

    args = parser.parse_args()

    Path(args.results_dir).mkdir(parents=True, exist_ok=True)

    params_embedder = {"batch_size": 512,
                    "num_workers": 4,
                    "device": "cuda",
                    "engine": "tiffslide",
                    "patch_size": args.patch_size,
                    "patch_level": 0,
                    "patch_downsample": 1.0}

    slide_embedder = ptb_emb.SlideEmbedder(slide_dir=args.directory_HE_reg,
                                           coords_dir=args.directory_patches, **params_embedder)

    slides_IHC = os.listdir(args.directory_IHC_reg)

    CELL_CLASSES = ["Neoplastic", "Inflammatory", "Connective", "Dead", "Epithelial"]
    
    if args.tdx_to_remove is not None:
        df = pd.read_csv(args.tdx_to_remove)
        tdx_nb_to_remove = df["tdx_nb"].values
        
    for slide in os.listdir(args.directory_HE_reg):
        slide_embedder.slide_dir = args.directory_HE_reg
        slide_name = slide.split(f".{args.slide_ext}")[0]
        tdx_nb = slide_name.split("-B0")[0]
        if args.tdx_to_remove is not None and tdx_nb in tdx_nb_to_remove:
            continue
        # Extract the embeddings
        slide_embedder.extract_cell_based_embeddings(slide_name=slide_name, slide_ext=args.slide_ext,
                                                    cells_path=os.path.join(args.directory_cell_detection_reg,
                                                                            f"{slide_name}.geojson"),
                                                    cell_classes=CELL_CLASSES,
                                                    with_offset=True)
        matching_slide_IHC = slides_IHC[np.where([tdx_nb in s for s in slides_IHC])[0].item()].split(f".{args.slide_ext}")[0]
        coords, attrs = read_h5_file(os.path.join(args.directory_patches, f"{matching_slide_IHC}.h5"), 'coords')
        slide_embedder.slide_dir = args.directory_IHC_reg
        slide_embedder.extract_stain_based_embeddings(matching_slide_IHC, slide_ext=args.slide_ext, coords=coords,
                                                      conv_matrix_name="HD_custom")
        slide_embedder.stain_based_embeddings[slide_name] = slide_embedder.stain_based_embeddings.pop(matching_slide_IHC)
        slide_embedder.save_embeddings(args.results_dir, merge=True, flush_memory=True)
        slide_embedder.save_embeddings_names(args.results_dir, merge=True, flush_memory=True)