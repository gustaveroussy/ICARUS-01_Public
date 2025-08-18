# Machine-Learning-based Analysis of Tumor Microenvironment and Target Antigen Spatial Distribution

This repository provides scripts and documentation for reproducing the AI-based analysis of the tumor microenvironment and target antigen spatial distribution conducted in the Icarus BREAST01 and Icarus LUNG01 clinical trials.

## Installation

To use this code, we recommend setting up a Conda environment with the required dependencies. The code was tested with Python 3.10 To install the dependencies manually, follow these steps:

1. **Create a Conda environment with Python 3.10:**

    ```bash
    conda create -n icarus python=3.10
    conda activate icarus
    ```

2. **Install the main dependencies using pip:**

    ```bash
    pip install valis-wsi==1.0.4
    pip install "prismtoolbox[emb,seg]==0.1.1"
    pip install opencv-contrib-python jupyter seaborn statannotations
    ```

To ensure reproducibility, we also provide an `environment.yml` and a `requirements.txt` file to replicate the exact environment used during our experiments. You can set up the environment using the following commands:

```bash
conda env create -f environment.yml
conda activate icarus
pip install -r requirements.txt --no-deps
```

For nuclei segmentation and classification on H&E slides, we used the [CellVit](https://github.com/TIO-IKIM/CellViT) model. Please follow the instructions in their repository to retrieve nuclei contours and classification as GeoJSON files.

## Usage

The Icarus project consists of several Python, Groovy, and R scripts organized as follows:

### Python Scripts (`python_scripts` folder)

1. `1_script_contouring.py`: Contours tissue in slides using [PrismToolbox](https://github.com/gustaveroussy/PrismToolBox) and saves the results as pickle files.
2. `2_script_patching.py`: Extracts patches from tissue contours and saves them as h5 files.
3. `3_script_embedding.py`: Creates embeddings from patches using pre-trained models and saves them as pt files.
4. `4_script_nuclei.py`: Performs nuclei segmentation using the [SOP](https://github.com/loic-lb/Unsupervised-Nuclei-Segmentation-using-Spatial-Organization-Priors) model and saves results as GeoJSON files.
5. `5_script_prepare_materials_for_reg.py`: Prepares data for the registration task.
6. `6_script_perform_registration.py`: Registers slides using the [VALIS](https://valis.readthedocs.io/en/latest/index.html) registration pipeline, saving registration models as pickle files.
7. `7_script_warp_coords.py`: Warps patch coordinates from one modality to another using the trained registration models.
8. `8_script_embedding_reg.py`: Computes embeddings based on nuclei segmentation from CellViT and staining features from registered patches.

### Groovy Scripts (`qupath_scripts` folder)

- **Nuclei:** Imports nuclei segmentation into QuPath, estimates cell membranes, and measures stain-based features.
- **Clustering:** Imports ROI from Python scripts into QuPath and saves them as OME-TIFF files.

### Notebooks (`notebooks` folder)

Notebooks are provided for processing data produced by the Python scripts and performing the analysis.

### R Scripts (`R_scripts` folder)

R scripts are used for generating figures and conducting statistical analyses.

## Reproducibility

To reproduce the Icarus project results, follow these steps:

1. **DAISY Clinical Trial-Like Analysis:**

    Run the following scripts to perform tissue contouring, patch extraction, and patch embedding:

    ```bash
    python python_scripts/1_script_contouring.py --slide_directory path_to_IHC_slides --slide_type IHC --annotations_directory path_to_annotations --results_directory path_to_results

    python python_scripts/2_script_patching.py --slide_directory path_to_IHC_slides --slide_type IHC --contours_directory path_to_IHC_contours --patch_size 80 --overlap 0 --results_directory path_to_results

    python python_scripts/3_script_embedding.py --slide_directory path_to_IHC_slides --patch_directory path_to_80x80_patches --arch_name clam_or_conch --results_directory path_to_results 
    ```

2. **Intracellular Stain Distribution-Based Analysis:**

    Run the `1_script_contouring.py` without annotations and `2_script_patching.py` with patch size 512 and 20 overlap, then perform nuclei segmentation using the SOP model:

    ```bash
    python python_scripts/4_script_nuclei.py --slide_directory path_to_IHC_slides --patch_directory path_to_512x512_patches --pretrained_weights path_to_pretrained_SOP_model --results_directory path_to_results
    ```
    You can retrieve the pretrained weights for the SOP model [here](https://github.com/loic-lb/Unsupervised-Nuclei-Segmentation-using-Spatial-Organization-Priors/tree/DAISY?tab=readme-ov-file). If performing BREAST analysis, add the `--clip_custom` tag to clamp staining values. After running this script, import the GeoJSON files into QuPath and use the Groovy scripts in the `qupath_scripts/nuclei` folder to generate cell membranes and measure stain-based features.

3. **Registration of the H&E and IHC Slides:**

    To analyze jointly the H&E and IHC slides through the registration of both modalities, first contour the tissue in each slide type using the `1_script_contouring.py` script (omit annotations for this step). Run the CellViT model on the H&E slides to retrieve nuclei contours and classifications. Follow CellViT’s repository instructions for details.

    Next, prepare data for the registration task:

    ```bash
    python python_scripts/5_script_prepare_materials_for_reg.py --directory_IHC path_to_IHC_slides --directory_HE path_to_H&E_slides --contours_IHC_dir path_to_IHC_contours --contours_HE_dir path_to_HE_contours --cellvit_detection_dir path_to_cellvit_seg --contours_to_keep path_to_contours_to_keep --results_directory path_to_results
    ```
    This script prepares data for registration. This script outputs:
    * A bounding box for each slide delineating the ROI that will be considered for registration (according to the contours indicated in `contours_to_keep`)
    * The nuclei segmentation from CellViT restricted to the contours present in the ROI defined by the bounding box, with their coordinates relative to the coordinates of the origin of the bounding box
    * The contours restricted to the ROI for the H&E and IHC slides
    
    You may then import the resulting ROIs into QuPath and save them as OME-TIFF slides using the Groovy scripts present in the `qupath_scripts/clustering` folder for both H&E and IHC slides.

    Follow with patch extraction on the OME-TIFF H&E slides using the script `2_script_patching.py` and a patch size of 128 with no overlap.
    
    Then, perform registration using the script `6_script_perform_registration.py` using the OME-TIFF slides and the contours restricted to the ROIs produced by the script `5_script_prepare_materials_for_reg.py`. This script will output registration models that can be used to warp patches from one modality to another.

    You can then warp patch coordinates you extracted from the H&E OME-TIFF slides to the IHC OME-TIFF slides using the script `7_script_warp_coords.py` (using the contours restricted to the ROIs from the IHC ROIs as target).

    Finally, compute embeddings based on the nuclei segmentation from CellViT processed by the script `5_script_prepare_materials_for_reg.py`, as well as staining features from the registered patches on IHC OME-TIFF slides, using the script `8_script_embedding_reg.py`.

    For ICARUS LUNG, the registration was also used to provide an interpretation of clustering results after DAISY-like analysis. You may reproduce this analysis by employing the same steps above, but, this time, using the  `7_script_warp_coords.py` script with the IHC slides as the source and the H&E slides as the target (you can keep the same registation models used for H&E to IHC warping). You also need to use the tag `--bbox_origin_dir` to offset the coordinates of the patches computed on the whole slide (of size 80) to the coordinates of the ROI of the IHC slide. You can then use the script `8_script_embedding_reg.py` to compute interpretable embeddings (be careful about the patch size that should be set to 80 here and not default 128).

For all these analyses, you may use the provided notebooks for data processing and the R scripts for figures and statistical analysis.

## Contact

For questions or inquiries, please contact the principal investigators. If you encounter any issues with the code, please open an issue in this repository.
