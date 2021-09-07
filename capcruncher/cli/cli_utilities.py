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




