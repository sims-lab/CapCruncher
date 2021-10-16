import click
from capcruncher.cli import UnsortedGroup

@click.group()
def cli():
    """Contains miscellaneous functions"""

@cli.command()
@click.argument('gtf')
@click.option('-o', '--output', help='Output file name')
def gtf_to_bed12(gtf: str, output:str):
    from pybedtools import BedTool
    from capcruncher.utils import gtf_line_to_bed12_line

    bt_gtf = BedTool(gtf)
    df_gtf = bt_gtf.to_dataframe()
    df_gtf['geneid'] = df_gtf['attributes'].str.extract(r'gene_id\s?\"(.*?)\";.*')
    df_gtf = df_gtf.query('feature.isin(["5UTR", "3UTR", "exon"])')
    df_gtf = df_gtf.loc[lambda df: df['seqname'].str.contains(r'^chr[xXYy]?[1-9]?[0-9]?$')]

    with open(output, 'w') as w:
        for gene, df in df_gtf.sort_values(['seqname', 'start']).groupby('geneid'):
            w.write(gtf_line_to_bed12_line(df) + '\n')


@cli.command()
@click.argument('bincode_file')
@click.argument('output_file')
def bincode_to_json(bincode_file, output_file):

    from libcapcruncher.fastq_deduplication import bincode_to_dict
    import ujson
    import xopen

    bc_dict = bincode_to_dict(bincode_file)

    with xopen.xopen(output_file, 'w') as w:
        ujson.dump(bc_dict, w)




