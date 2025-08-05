import argparse
import numpy as np
import pandas as pd
import re
import functools
print = functools.partial(print, flush=True)

DataFrame = pd.core.frame.DataFrame

def explode_df(df: DataFrame, cols: list, sep=',', fill_value: str='', preserve_index: bool=False) -> DataFrame:
    """
    Expand dataframe entries of the columns specified in l_cols and for which there are multiple values.

    Parameters
    ---------
    df: DataFrame
        Input dataframe on which expansion is performed
    cols: list
        List of columns where expansion is required
    sep: char
        Character separating the multiple values
        Default : ','
    fill_value: bool
        Entry in exanpded dataframe for empty lists
        Default : ''
    preserve_index: bool
        Whether original index should be preserved or not. If set to True, the index of the expanded DataFrame
        will be redundant.
        Default : False

    Returns
    -------
    df_expanded: DataFrame
        Returns  dataframe where entries of any of the columns in l_cols with multiple values have been expanded.
    """
    # transform comma-separated to list
    df = df.assign(**{col:df[col].str.split(sep) for col in cols}).copy()
    if (cols is not None and len(cols) > 0 and not isinstance(cols, (list, tuple, np.ndarray, pd.Series))):
        cols = [cols]
    # calculate lengths of lists
    lens = df[cols[0]].str.len()
    # format NaN to [NaN] and strip unwanted characters
    for col in cols:
        df.loc[df[col].isnull(), col] = df.loc[df[col].isnull(), col].apply(lambda x: [np.nan])
        df.loc[lens > 1, col] = df.loc[lens > 1, col].apply(lambda x: [y.strip() for y in x])
    # all columns except `cols`
    idx_cols = df.columns.difference(cols)
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    df_expanded = (pd.DataFrame({col:np.repeat(df[col].values, lens) for col in idx_cols},
                index=idx).assign(**{col:np.concatenate(df.loc[lens>0, col].values) for col in cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        df_expanded = (df_expanded.append(df.loc[lens==0, idx_cols], sort=False).fillna(fill_value))
    # revert the original index order
    df_expanded = df_expanded.sort_index()
    # reset index if requested
    if not preserve_index:
        df_expanded = df_expanded.reset_index(drop=True)
    return df_expanded


def print_status_changes(df_table):
    old_ok = ~df_table["IRODS_Status_Old"].isnull()
    new_ok = ~df_table["IRODS_Status"].isnull()
    total = df_table.shape[0]
    changed_to_ok = ~old_ok & new_ok
    changed_to_nok = old_ok & ~new_ok
    print("%d/%d samples changed from Not OK to OK:" % (sum(changed_to_ok), total))
    if sum(changed_to_ok)>0:
        sample_ids_changed_to_ok = df_table.loc[changed_to_ok, "Sample_Id"].tolist()
        print("\t" + "\n\t".join(sample_ids_changed_to_ok))
    print("%d/%d samples changed from OK to Not OK" % (sum(changed_to_nok), total))
    if sum(changed_to_nok)>0:
        sample_ids_changed_to_nok = df_table.loc[changed_to_nok, "Sample_Id"].tolist()
        print("\t" + "\n\t".join(sample_ids_changed_to_nok))
    return sum(changed_to_ok), sum(changed_to_nok)


def main(args):
    df_table = pd.read_table(args.table)
    if "IRODS_Status" not in df_table:
        raise ValueError("cannot check IRODS status as 'IRODS_Status' column is absent from %s" % args.table)
    else:
        df_table = df_table.rename(columns={"IRODS_Status": "IRODS_Status_Old"})

    df_irods = pd.read_table(args.irods, header=None)
    df_irods.columns = ["Dataset"]
    mask_keep = df_irods["Dataset"].apply(lambda x: x.startswith("collection:"))
    df_irods = df_irods.loc[mask_keep].copy()
    df_irods["Dataset"] = df_irods["Dataset"].apply(lambda x: re.sub("collection: ", "", x))
    df_irods["IRODS_Status"] = "OK"

    # explode possible comma-separated values of FASTQ_1/2
    df_table = explode_df(df=df_table, cols=["FASTQ_1", "FASTQ_2", "FASTQ_1_Name", "FASTQ_2_Name"], sep=",")

    col_x = "FASTQ_1"
    col_y = "Dataset"
    df_table[col_y] = df_table[col_x].str.extract(r'^(.+)(?=\/archive)')
    df_table_1 = df_table.merge(df_irods, how="left", on="Dataset")

    col_x = "FASTQ_2"
    df_table[col_y] = df_table[col_x].str.extract(r'^(.+)(?=\/archive)')
    df_table_2 = df_table.merge(df_irods, how="left", on="Dataset")

    n_changed_to_ok_1, n_changed_to_nok_1 = print_status_changes(df_table_1)
    n_changed_to_ok_2, n_changed_to_nok_2 = print_status_changes(df_table_2)

    if n_changed_to_ok_1 + n_changed_to_nok_1 + n_changed_to_ok_2 + n_changed_to_nok_2 > 0:
        raise ValueError("some files changed status after automatic IRODS check. Please review your " + \
                         "the values of 'IRODS_Status' in your config file 'samples.tsv' or disable automatic " + \
                         "IRODS check in the config.yaml file.")
    else:
        open(args.output, 'w').close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Checks whether all files dataset are still
                                     available in IRODS or not""")
    parser.add_argument("--table", type=str, help="Path to table of samples")
    parser.add_argument("--irods", type=str, help="Path to output of imeta qu -C command.")
    parser.add_argument("--output", type=str, help="Path to output table with status.")

    args = parser.parse_args()
    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
