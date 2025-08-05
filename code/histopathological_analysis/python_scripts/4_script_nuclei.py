import os
import argparse
from pathlib import Path

import prismtoolbox.nucleiseg as ptb_nucleiseg



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract nuclei from patches")
    parser.add_argument("--slide_directory", type=str, help="Path to the directory containing the files")
    parser.add_argument("--patch_directory", type=str, help="Path to the directory containing the patches")
    parser.add_argument("--pretrained_weights", type=str, 
                        help="Path to the pretrained weights for the nuclei segmentation model")
    parser.add_argument("--slide_ext", type=str, default="svs", help="Extension of the slide files")
    parser.add_argument("--clip_custom", action="store_true", help="Whether to clip the patches with custom values")
    parser.add_argument("--clip_min_value", type=float, default=0.3, help="Minimum value for clipping the patches")
    parser.add_argument("--clip_max_value", type=float, default=1.0, help="Maximum value for clipping the patches")
    parser.add_argument("--results_directory", type=str, default=None, help="Path to the directory where the results will be saved")
    args = parser.parse_args()
    
    Path(args.results_directory).mkdir(parents=True, exist_ok=True)

    # Use the following parameters to extract nuclei from the patches
    params_nucleiseg = {"model_name": "sop",
                        "threshold_overlap": 0.5,
                        "pretrained_weights": args.pretrained_weights,
                        "batch_size": 10,
                        "num_workers": 4,
                        "device": "cuda"}
    
    if args.clip_custom:
        transforms_dict = {"totensor": {}, 
                           "clip_custom": {"min_value": args.clip_min_value, "max_value": args.clip_max_value},
                           "normalize": {"mean": (0.5, 0.5, 0.5), "std": (0.5, 0.5, 0.5)}}
    else:
        transforms_dict = {"totensor": {}, 
                           "normalize": {"mean": (0.5, 0.5, 0.5), "std": (0.5, 0.5, 0.5)}}

    slide_embedder = ptb_nucleiseg.NucleiSegmenter(slide_dir=args.slide_directory,
                                                   coords_dir=args.patch_directory, 
                                                   transforms_dict=transforms_dict,
                                                   **params_nucleiseg)

    for slide in os.listdir(args.slide_directory):
        slide_name = slide.split(f".{args.slide_ext}")[0]
        print(slide_name)
        if f"{slide_name}.geojson" in os.listdir(args.results_directory):
            continue
        # Extract the nuclei
        slide_embedder.segment_nuclei(slide_name, slide_ext=args.slide_ext)
        # Save the nuclei
        slide_embedder.save_nuclei(args.results_directory, slide_ext=args.slide_ext)
        break