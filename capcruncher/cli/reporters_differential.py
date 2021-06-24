import os
import pandas as pd
import itertools
import numpy as np


def get_chromosome_from_name(df: pd.DataFrame, name: str):
    chrom = (df.query(f'name == "{name}"')
               ['chrom']
               .iloc[0])
    return chrom
    

def differential(union_bedgraph: os.PathLike,
                              capture_name: str,
                              capture_viewpoints: os.PathLike,
                              output_prefix: os.PathLike = 'differential',
                              design_matrix: os.PathLike = None,
                              grouping_col: str = 'condition',
                              threshold_count: float = 20,
                              threshold_q: float = 0.05,
                              threshold_mean: float = 0):
    
    """
    Identifies differential interactions between conditions.

    Parses a union bedgraph containg reporter counts from at least two conditions with
    two or more replicates for a single capture probe and outputs differential interaction
    results. Following filtering to ensure that the number of interactions is above the required 
    threshold (--threshold_count), diffxpy is used to run a wald test after 
    fitting a negative binomial model to the interaction counts.The options to filter
    results can be filtered by a minimum mean value (threshold_mean) and/or 
    maximum q-value (threshold-q) are also provided.

    Notes:
        
     Currently both the capture oligos and the name of the probe being analysed must 
     be provided in order to correctly extract cis interactions.
 
     If a N_SAMPLE * METADATA design matrix has not been supplied, the script 
     assumes that the standard replicate naming structure has been followed 
     i.e. SAMPLE_CONDITION_REPLICATE_(1|2).fastq.gz.    

    \f
    Args:
        union_bedgraph (os.PathLike): Union bedgraph containg all samples to be compared.
        capture_name (str): Name of capture probe. MUST match one probe within the supplied oligos.
        capture_viewpoints (os.PathLike): Capture oligos used for the analysis.
        output_prefix (os.PathLike, optional): Output prefix for differntial interactions. Defaults to 'differential'.
        design_matrix (os.PathLike, optional): Design matrix to use for grouping samples. (N_SAMPLES * METADATA). Defaults to None.
        grouping_col (str, optional): Column to use for grouping. Defaults to 'condition'.
        threshold_count (float, optional): Minimum number of reported interactions required. Defaults to 20.
        threshold_q (float, optional): Maximum q-value for output. Defaults to 0.05.
        threshold_mean (float, optional): Minimum mean value for output. Defaults to 0.
    """    

    import diffxpy.api as de
    
    df_bdg = pd.read_csv(union_bedgraph, sep='\t')
    df_viewpoints = pd.read_csv(capture_viewpoints, sep='\t', names=['chrom', 'start', 'end', 'name'])

    #  If design matrix present then use it. Else will assume that the standard format has been followed:
    #  i.e. NAME_TREATMENT_REPLICATE
    if design_matrix:
        df_design = pd.read_csv(design_matrix, sep='\t')
    else:
        col_dict = {col: '_'.join(col.split('_')[:-1]) for col in df_bdg.columns[3:]}
        df_design = pd.Series(col_dict).to_frame(grouping_col)
    

    # Only cis interactions
    capture_chrom = get_chromosome_from_name(df_viewpoints, name=capture_name)
    df_bdg_counts = df_bdg.query(f'chrom == "{capture_chrom}"')

    # Only counts
    df_bdg_counts = df_bdg_counts.iloc[:, 3:]

    # Only with number of interactions > threshold per group in at least 2 replicates
    df_bdg_counts = (df_bdg_counts.groupby(df_design[grouping_col], axis=1)
                                  .apply(lambda df: df[(df >= threshold_count).sum(axis=1) >= 2])
                                  .fillna(0.0))
    
    # Run differential testing
    count_data = df_bdg_counts.transpose().values.astype(np.float64)
    fragment_names = df_bdg_counts.index.values

    tests = de.test.pairwise(count_data, 
                             grouping=grouping_col, 
                             sample_description=df_design,
                             gene_names=fragment_names,
                             test='wald',
                             lazy=False, 
                             backend='numpy')
       
    # Go through all of the pairwise tests
    for g1, g2 in itertools.combinations(tests.groups, 2):
        df_comparison = tests.summary_pairs(groups0=[g1,], 
                                            groups1=[g2,],
                                            qval_thres=threshold_q, 
                                            mean_thres=threshold_mean)

        # Extract the fragment coords
        df_coords = df_bdg.loc[df_comparison['gene'], ['chrom', 'start', 'end']]
        # Merge with test results
        df_comparison = df_comparison.merge(df_coords, left_on='gene', right_index=True)
        # Output to tsv
        (df_comparison.drop(columns='gene')
                      [['chrom', 'start', 'end', 'mean', 'log2fc', 'pval', 'qval']]
                      .to_csv(f'{output_prefix}_{g1}_vs_{g2}.tsv', sep='\t', index=False))



# if __name__ == '__main__':
#     interactions_differential()
        




    
