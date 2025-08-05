import os
import argparse
import logging
import sys
from pathlib import Path

import prismtoolbox as ptb

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract patches from slides with contours already extracted")
    parser.add_argument("--slide_directory", type=str, help="Path to the directory containing the slides")
    parser.add_argument("--slide_type", type=str, choices=["HE", "IHC"], help="Slide type between 'HE' and 'IHC'")
    parser.add_argument("--contours_directory", type=str, help="Path to the directory containing the contours")
    parser.add_argument("--engine", type=str, default="openslide",choices=["openslide", "tiffslide"],
                        help="Engine to use for reading the slides")
    parser.add_argument("--arch_name", type=int, default=128, help="Size of the patches to extract (in pixels)")
    parser.add_argument("--overlap", type=int, default=0, help="Overlap between patches (in micrometers)")
    parser.add_argument("--results_directory", type=str, default=None, help="Path to the directory where the results will be saved")
    args = parser.parse_args()
    
    params_patches = {"patch_size": args.patch_size, "patch_level": 0, "overlap": args.overlap, 
                      "mode": "contours", "units": ("px", "micro"), "contours_mode": "four_pt", "rgb_threshs": (2, 240)}
    
    # Use the following parameters to stitch the extracted patches
    params_stitch_WSI = {"vis_level": 2, "draw_grid": False}
    # Set the suffix to add to the directory names
    suffix = "_annotations" if "annotations" in args.contours_directory else ""

    # Path to the directory where the patches will be saved as h5 files
    directory_patches = os.path.join(args.results_directory,
                                     f"patches_{params_patches['patch_size']}_overlap"
                                     f"_{params_patches['overlap']}{suffix}_{args.slide_type}")
    # Path to the directory where the patches will be saved as geojson
    directory_patches_qupath = os.path.join(args.results_directory,
                                            f"patches_{params_patches['patch_size']}_overlap"
                                            f"_{params_patches['overlap']}_qupath{suffix}_{args.slide_type}")
    # Path to the directory where the stitched images will be saved as jpg files
    directory_stitch = os.path.join(args.results_directory,
                                    f"stitch_images_{params_patches['patch_size']}_overlap"
                                    f"_{params_patches['overlap']}{suffix}_{args.slide_type}")
    
    Path(directory_patches).mkdir(parents=True, exist_ok=True)
    #Path(directory_patches_qupath).mkdir(parents=True, exist_ok=True)
    Path(directory_stitch).mkdir(parents=True, exist_ok=True)

    # Iterate over the files in the directory
    for file_name in os.listdir(args.slide_directory):
        # Load the image
        WSI_object = ptb.WSI(os.path.join(args.slide_directory, file_name), engine=args.engine)
        print(f"Processing {WSI_object.slide_name}...")

        if f"{WSI_object.slide_name}.h5" in os.listdir(directory_patches):
            continue
        
        # Load the contours for the image
        WSI_object.load_tissue_contours(os.path.join(args.contours_directory, f"{WSI_object.slide_name}.pkl"))

        if len(WSI_object.tissue_contours) == 0:
            WSI_object.set_roi()
            params_patches["mode"] = "roi"
        else:
            params_patches["mode"] = "contours"

        # Extract patches from the contours
        WSI_object.extract_patches(**params_patches)
        # Save the extracted patches as h5 files
        WSI_object.save_patches(directory_patches, file_format="h5")
        # Save the extracted patches as geojson files
        #WSI_object.save_patches(directory_patches_qupath, file_format="geojson")

        # Stitch the extracted patches
        if params_stitch_WSI["vis_level"] >= len(WSI_object.level_dimensions):
            continue
        img = WSI_object.stitch(**params_stitch_WSI)

        img.save(os.path.join(directory_stitch, f"{WSI_object.slide_name}.jpg"))