import os
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import prismtoolbox as ptb
from shapely import MultiPolygon, Polygon
from tqdm import tqdm

COLOR_DICT = {
    1: [255, 0, 0],
    2: [34, 221, 77],
    3: [35, 92, 236],
    4: [254, 255, 0],
    5: [255, 159, 68],
}

TYPE_NUCLEI_DICT = {
    1: "Neoplastic",
    2: "Inflammatory",
    3: "Connective",
    4: "Dead",
    5: "Epithelial",
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare materials for registration")
    parser.add_argument("--result_dir", type=str, 
                        help="Path to the directory where the results will be saved")
    parser.add_argument("--directory_IHC", type=str, help="Path to the directory containing the IHC slides")
    parser.add_argument("--directory_HE", type=str, help="Path to the directory containing the H&E slides")
    parser.add_argument("--contours_IHC_dir", type=str, help="Path to the directory containing the contours extracted from IHC slides")
    parser.add_argument("--contours_HE_dir", type=str, help="Path to the directory containing the contours extracted from H&E slides")
    parser.add_argument("--cellvit_detection_dir", type=str, 
                        help="Path to the directory containing the cell detection results from CellVit")
    parser.add_argument("--contours_to_keep", type=str, help="Path to the csv file containing the contours to keep for each slide")
    parser.add_argument("--annotation_directory", type=str, default=None, help="Path to the directory containing the pathologist annotations")
    
    args = parser.parse_args()
    
    slides_IHC = os.listdir(args.directory_IHC)
    slides_HE = os.listdir(args.directory_HE)
    
    contours_df = pd.read_csv(args.contours_to_keep)
    
    cell_detection_save_path = os.path.join(args.result_dir, "cell_detection_reg")
    Path(cell_detection_save_path).mkdir(parents=True, exist_ok=True)
    
    for slide_type in ["HE", "IHC"]:
        bbox_save_path = os.path.join(args.result_dir, f"bbox_{slide_type}_reg")
        contours_save_path = os.path.join(args.result_dir, f"contours_{slide_type}_reg")
        Path(bbox_save_path).mkdir(parents=True, exist_ok=True)
        Path(contours_save_path).mkdir(parents=True, exist_ok=True)
        
    print("Processing the slides and cell detections for registration...")
    for slide_IHC in tqdm(slides_IHC):
        tdx_nb = slide_IHC.split("-B0")[0]
        slide_HE = slides_HE[np.where([slide_IHC.split("-B0")[0] in s for s in slides_HE])[0].item()]
        slide_IHC_path = os.path.join(args.directory_IHC, slide_IHC)
        slide_HE_path = os.path.join(args.directory_HE, slide_HE)
        contours_IHC_path = os.path.join(args.contours_IHC_dir, f"{slide_IHC.split('.svs')[0]}.pkl")
        contours_HE_path = os.path.join(args.contours_HE_dir, f"{slide_HE.split('.svs')[0]}.pkl")
        contours_to_keep_df = contours_df[contours_df["TdxNumber"] == tdx_nb]
        for slide_type, slide_path, contours_path in zip(["HE", "IHC"], [slide_HE_path, slide_IHC_path],
                                                            [contours_HE_path, contours_IHC_path]):
            # Load the slide and the contours
            WSI_object = ptb.WSI(slide_path, engine="openslide")
            WSI_object.load_tissue_contours(contours_path)
            if args.annotation_directory is not None:
                WSI_object.apply_pathologist_annotations(os.path.join(args.annotation_directory, f"{WSI_object.slide_name}.geojson"))
            contours_to_keep = list(map(int, contours_to_keep_df[f"contours_{slide_type}"].values[0].split(" ")))
            # Keep only the contours to keep
            WSI_object.tissue_contours = [WSI_object.tissue_contours[contours_idx] for contours_idx in contours_to_keep]
            # Extract the bounding box of the contours
            xmin, ymin, xmax, ymax = ptb.utils.vis_utils.bbox_from_contours(WSI_object.tissue_contours)
            # Create a polygon from the bounding box
            polygons = [Polygon([
                        (xmin, ymin),
                        (xmax, ymin),
                        (xmax, ymax),
                        (xmin, ymax),
                    ])]
            polygons = MultiPolygon(polygons)
            bbox_save_path = os.path.join(args.result_dir, f"bbox_{slide_type}_reg")
            # Save the bounding box as geojson for use in QuPath
            ptb.utils.qupath_utils.export_polygons_to_qupath(polygons, os.path.join(bbox_save_path, f"{WSI_object.slide_name}.geojson"),
                                                             "annotation", offset=WSI_object.offset)
            if slide_type == "HE":
                # Load the cell detections from CellVit
                cells = ptb.utils.data_utils.load_obj_with_json(os.path.join(args.cellvit_detection_dir, WSI_object.slide_name,
                                                                            "cell_detection/cells.json"))
                cells_list = cells["cells"]
                # Create a dataframe from the cell detections
                cell_segmentation_df = pd.DataFrame(cells_list)
                detected_types = sorted(cell_segmentation_df.type.unique())
                # Create a polygon from the tissue contours
                contours_polygons = ptb.utils.qupath_utils.contoursToPolygons(WSI_object.tissue_contours)
                # Set the offset
                offset = (-xmin, -ymin)
                for cell_type in detected_types:
                    # Keep only the cells of the current type
                    cells = cell_segmentation_df[cell_segmentation_df["type"] == cell_type]
                    contours = cells["contour"].to_list()
                    # Create a polygon from the cell contours
                    cell_polygons = ptb.utils.qupath_utils.contoursToPolygons([np.array(cnt) for cnt in contours], make_valid=True)
                    # Get the intersection between the cell polygons and the tissue contours
                    polygons_in_contours = ptb.utils.qupath_utils.intersectionPolygons(cell_polygons, contours_polygons)
                    ptb.utils.qupath_utils.export_polygons_to_qupath(polygons_in_contours, 
                                                                    os.path.join(cell_detection_save_path,
                                                                                 f"{WSI_object.slide_name}.geojson"), 
                                                                    "detection",
                                                                    label=TYPE_NUCLEI_DICT[cell_type],
                                                                    color=COLOR_DICT[cell_type],
                                                                    append_to_existing_file=True,
                                                                    offset=offset)
          
            # Crop the tissue contours to the bounding box
            WSI_object.tissue_contours = [(np.array(contour) - (xmin, ymin)).astype("int64") for contour in WSI_object.tissue_contours]
            # Save the tissue contours
            contours_save_path = os.path.join(args.result_dir, f"contours_{slide_type}_reg")
            WSI_object.save_tissue_contours(contours_save_path)