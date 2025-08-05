# -*- coding: utf-8 -*-
"""
@created: Aug 21 2024
@modified: Aug 21 2024
@author: Yoann Pradat

    Institut Gustave Roussy
    U981 & Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Prepare a table summarizing all the data available and analyzed in the DAISY project.
"""

import argparse
import os
import numpy as np
import pandas as pd
import random

# functions ============================================================================================================

def main(args):
    # load table listing all patients
    df_ids_clear = pd.read_table(args.ids_clear)

    # anonymize breast and lung cohorts separately
    mask_b = df_ids_clear["Cancer_Type"]=="Breast Cancer"
    mask_l = df_ids_clear["Cancer_Type"]=="Lung Cancer"
    n_ids_b = sum(mask_b)
    n_ids_l = sum(mask_l)
    df_ids_anonym_b = df_ids_clear.loc[mask_b, ["Subject_Id"]].copy()
    df_ids_anonym_l = df_ids_clear.loc[mask_l, ["Subject_Id"]].copy()

    # randomly shuffle breast
    ids_b_shuffled = [i+1 for i in range(n_ids_b)]
    random.shuffle(ids_b_shuffled)
    df_ids_anonym_b["Subject_Id_Anonymous"] = ["ICARUS-B01-%03d" % x for x in ids_b_shuffled]

    # randomly shuffle lung
    ids_l_shuffled = [i+1 for i in range(n_ids_l)]
    random.shuffle(ids_l_shuffled)
    df_ids_anonym_l["Subject_Id_Anonymous"] = ["ICARUS-L01-%03d" % x for x in ids_l_shuffled]

    # assemble
    df_ids_anonym = pd.concat((df_ids_anonym_b, df_ids_anonym_l))

    # save
    df_ids_anonym.to_csv(args.ids_anonym, sep="\t", index=False)
    print("-file saved at %s" % args.ids_anonym)


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Prepare a table anonymizing patient ids.')
    parser.add_argument('--ids_clear', type=str, help='Path to table with all patient ids.',
                      default="../../data/cln/curated/cln_icarus_curated_with_qc.tsv")
    parser.add_argument('--ids_anonym', type=str, help='Path to output anonymisation table.',
                        default="../../results/tables_paper/table_anonymisation.tsv")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
