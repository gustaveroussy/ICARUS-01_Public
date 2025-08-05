import sys
import os
import logging
import argparse
import prismtoolbox.wsiemb as ptb_emb
from pathlib import Path

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract embeddings from slides with patches already extracted")
    parser.add_argument("--slide_directory", type=str, help="Path to the directory containing the slides")
    parser.add_argument("--engine", type=str, default="openslide", choices=["openslide", "tiffslide"],
                        help="Engine to use for reading the slides")
    parser.add_argument("--patch_directory", type=str, help="Path to the directory containing the patches")
    parser.add_argument("--patch_size", type=int, default=80, help="Size of the patches to extract (in pixels)")
    parser.add_argument("--slide_ext", type=str, default="svs", help="Extension of the slide files")
    parser.add_argument("--arch_name", type=str, default="clam", choices=["clam", "conch"],
                        help="Name of the architecture to use for the embeddings")
    parser.add_argument("--results_directory", type=str, default=None, help="Path to the directory where the results will be saved")
    args = parser.parse_args()

    # Use the following parameters to extract embeddings from the patches
    if args.arch_name == "clam":
        mean = [0.485, 0.456, 0.406]
        std = [0.229, 0.224, 0.225]
        pretrained_weights = "IMAGENET1K_V1"
        need_login=False
    elif args.arch_name == "conch":
        mean = [0.48145466, 0.4578275, 0.40821073]
        std = [0.26862954, 0.26130258, 0.27577711]
        pretrained_weights = None
        need_login=True
    params_embedder = {"arch_name": args.arch_name,
                       "pretrained_weights": pretrained_weights,
                       "transforms_dict": {"totensor": {}, "normalize": {"mean": mean, "std": std}},
                       "batch_size": 512,
                       "num_workers": 4,
                       "device": "cuda",
                       "engine": args.engine,
                       "need_login": need_login}

    # Path to the directory where the embeddings will be saved as pt files
    directory_embeddings = os.path.join(args.results_directory, f"./feats_{args.patch_size}_{args.arch_name}")
    Path(directory_embeddings).mkdir(parents=True, exist_ok=True)

    slide_embedder = ptb_emb.SlideEmbedder(slide_dir=args.slide_directory,
                                           coords_dir=args.patch_directory,
                                           **params_embedder)

    for slide in os.listdir(args.slide_directory):
        slide_name = slide.split(f".{args.slide_ext}")[0]
        if os.path.exists(os.path.join(directory_embeddings, f"{slide_name}.pt")):
            continue
        # Extract the embeddings
        slide_embedder.extract_model_based_embeddings(slide_name, slide_ext=args.slide_ext)
        slide_embedder.save_embeddings(directory_embeddings, flush_memory=True, merge=True)

