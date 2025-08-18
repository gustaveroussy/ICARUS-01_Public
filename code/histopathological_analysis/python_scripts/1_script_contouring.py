import os
import argparse
from pathlib import Path

import prismtoolbox as ptb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract tissue contours from slides")
    parser.add_argument("--slide_directory", type=str, help="Path to the directory containing the files")
    parser.add_argument("--slide_type", type=str, choices=["HE", "IHC"], help="Slide type between 'HE' and 'IHC'")
    parser.add_argument("--engine", type=str, choices=["openslide", "tiffslide"], default="openslide", help="Engine to use for reading the slides")
    parser.add_argument("--annotations_directory", type=str, default=None, 
                        help="Path to the directory containing the pathologist annotations")
    parser.add_argument("--results_directory", type=str, default=None, help="Path to the directory where the results will be saved")
    args = parser.parse_args()
    
    # Use the following parameters to detect tissue on the slides depending on the slide type
    if args.slide_type == "HE":
        params_detect_tissue = {"seg_level": 2, "window_avg": 30, "window_eng": 3, "thresh": 120, "area_min": 6e3}
    elif args.slide_type == "IHC":
        params_detect_tissue = {"seg_level": 2, "window_avg": 30, "window_eng": 5, "thresh": 90, "area_min": 6e3}
    else:
        raise ValueError("slide_type must be either 'HE' or 'IHC'")
    
        # Use the following parameters to visualize the extracted contours
    params_visualize_WSI = {"vis_level": 2, "number_contours": True, "view_slide_only": False, "line_thickness": 50}
    # Set the suffix to add to the directory names
    suffix = "_annotations" if args.annotations_directory is not None else ""

    # Path to the directory where the contours will be saved as pickle files
    directory_contours = os.path.join(args.results_directory, f"contours_{args.slide_type}{suffix}")
    # Path to the directory where the visualizations will be saved as jpg files
    directory_visualize = os.path.join(args.results_directory, f"contoured_images_{args.slide_type}{suffix}")

    Path(directory_contours).mkdir(parents=True, exist_ok=True)
    Path(directory_visualize).mkdir(parents=True, exist_ok=True)

    # Iterate over the files in the directory
    for file_name in os.listdir(args.slide_directory):
        # Load the image
        WSI_object = ptb.WSI(os.path.join(args.slide_directory, file_name), engine=args.engine)
        print(f"Processing {WSI_object.slide_name}...")

        if f"{WSI_object.slide_name}.pkl" in os.listdir(directory_contours):
            continue
        
        # Extract the contours from the image
        WSI_object.detect_tissue(**params_detect_tissue)
        # Apply pathologist annotations
        if args.annotations_directory is not None:
            WSI_object.apply_pathologist_annotations(os.path.join(args.annotations_directory, f"{WSI_object.slide_name}.geojson"))
        # Save extracted contours as pickle file
        WSI_object.save_tissue_contours(directory_contours)
        # Visualize the extracted contours on the tissue
        img = WSI_object.visualize(**params_visualize_WSI)
        img.save(os.path.join(directory_visualize, f"{WSI_object.slide_name}.jpg"))