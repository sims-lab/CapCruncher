import pandas as pd
import argparse


def main(rf_counts, rf_map, outfile=None):


    df_rf = pd.read_csv(rf_counts, header=None, names=['rf1', 'rf2', 'count'])
    df_rf = pd.read_csv(rf_map,
                    sep='\t',
                    header=None,
                    names=['chrom', 'start', 'end', 'name'])
    
    df_rf['mo_coord'] = (df_rf['chrom'].str.replace('chr', '') +
                     ':' +
                     df_rf['start'].astype(str) +
                     '-' +
                     df_rf['end'].astype(str))
    
    rf_name_to_coord_name = df_rf.set_index('name')['mo_coord'].to_dict()

    df_rf['rf1_coord'] = df_rf['rf1'].map(rf_name_to_coord_name)
    df_rf['rf2_coord'] = df_rf['rf2'].map(rf_name_to_coord_name)

    df_rf[['rf1_coord', 'rf2_coord', 'count']].to_csv(outfile, sep='\t', header=False, index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--rf_counts')
    parser.add_argument('-m', '--rf_map')
    parser.add_argument('-o', '--outfile', default='as_to_mo_converted.tsv.gz')
    args = parser.parse_args()

    main(**vars(args))


