import pandas as pd
import argparse


def main(rf_counts, rf_map, outfile=None):

    
    # Load rf counts
    df_rf_count = pd.read_csv(rf_counts,
                        header=None,
                        names=['rf1', 'rf2', 'count'],
                        )
    print('Loaded restriction fragment counts')
    
    # Load rf map
    df_rf_map = pd.read_csv(rf_map,
                    sep='\t',
                    header=None,
                    names=['chrom', 'start', 'end', 'name'])
    
    print('Loaded restriction fragment map')
    
    # Concatenate chrom start end for MO names
    df_rf_map['mo_coord'] = (df_rf_map['chrom'].str.replace('chr', '') +
                     ':' +
                     df_rf_map['start'].astype(str) +
                     '-' +
                     df_rf_map['end'].astype(str))
    
    # Generate a mapping
    rf_name_to_coord_name = df_rf_map.set_index('name')['mo_coord'].to_dict()
    print('Created mapping')

    # Map new name to counts
    df_rf_count['rf1_coord'] = df_rf_count['rf1'].map(rf_name_to_coord_name)
    df_rf_count['rf2_coord'] = df_rf_count['rf2'].map(rf_name_to_coord_name)

    # Output csv
    print('Outputting converted csv')
    df_rf_count[['rf1_coord', 'rf2_coord', 'count']].to_csv(outfile, sep='\t', header=False, index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--rf_counts')
    parser.add_argument('-m', '--rf_map')
    parser.add_argument('-o', '--outfile', default='as_to_mo_converted.tsv.gz')
    args = parser.parse_args()

    main(**vars(args))


