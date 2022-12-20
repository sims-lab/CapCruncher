import trackhub
import seaborn as sns
from capcruncher.utils import categorise_tracks
import pandas as pd
import os
import shutil

capcruncher_report = snakemake.input.report
bigwigs = [
    *snakemake.input.bigwigs_per_viewpoint,
    *snakemake.input.bigwigs_compared,
    *snakemake.input.bigwigs_summarised,
]
viewpoints = snakemake.input.viewpoints

hub_name = snakemake.config["hub"]["name"]

df_bigwigs = (
    pd.Series(bigwigs)
    .to_frame("fn")
    .assign(basename=lambda df: df["fn"].apply(os.path.basename))
)
attributes = df_bigwigs["basename"].str.extract(
    r"(?P<samplename>.*?)\.(?P<method>.*?)\.(?P<viewpoint>.*?)\.(?P<file_type>.*)"
)

df_bigwigs = (
    df_bigwigs.join(attributes)
    .assign(track_categories=lambda df: categorise_tracks(df["method"]))
    .sort_values(["samplename", "method", "viewpoint"])
)

# Create a hub
hub = trackhub.Hub(
    hub=hub_name,
    short_label=snakemake.config["hub"]["short"]
    if "short" in snakemake.config["hub"]
    else hub_name,
    long_label=snakemake.config["hub"]["long"]
    if "long" in snakemake.config["hub"]
    else hub_name,
    email=snakemake.config["hub"]["email"],
)


# Set up the genome
# This needs to be a custom genome if we are using a custom genome
if snakemake.config["genome"].get("custom"):
    genome = trackhub.Assembly(
        genome=snakemake.config["genome"]["name"],
        twobit_file=snakemake.config["genome"]["twobit"],
        organism=snakemake.config["genome"]["organism"],
        defaultPos=snakemake.config["hub"].get("default_position", "chr1:1000-2000"),
    )
    groups_file = trackhub.GroupsFile(
        [
            trackhub.GroupDefinition(
                name=hub_name, priority=1, default_is_closed=False
            ),
        ]
    )
    genome.add_groups(groups_file)

else:
    genome = trackhub.Genome(genome=snakemake.config["genome"]["name"])
    groups_file = None


# Create genomes file
genomes_file = trackhub.GenomesFile()

# Create trackdb
trackdb = trackhub.TrackDb()

# Add these to the hub
hub.add_genomes_file(genomes_file)
genome.add_trackdb(trackdb)
genomes_file.add_genome(genome)

# Extract groups for generating composite tracks
unique_samples = df_bigwigs["samplename"].unique()
unique_viewpoints = df_bigwigs["viewpoint"].unique()
unique_comparison_methods = df_bigwigs["method"].unique()

subgroup_vp = trackhub.SubGroupDefinition(
    name="viewpoint",
    label="Viewpoint",
    mapping={n.lower(): n for n in unique_viewpoints},
)
subgroup_sample = trackhub.SubGroupDefinition(
    name="samplename",
    label="Sample_Name",
    mapping={n.lower(): n for n in unique_samples},
)
subgroup_method = trackhub.SubGroupDefinition(
    name="summary_method",
    label="Summary_Method",
    mapping={n.split("-")[0]: n.split("-")[0] for n in unique_comparison_methods},
)

# Generate a color mapping based on sample names
colors = sns.color_palette("hls", len(unique_samples))
color_mapping = dict(zip(unique_samples, colors))

#######################
# # Add tracks to hub #
#######################

for category_name, df in df_bigwigs.groupby("track_categories"):

    composite = trackhub.CompositeTrack(
        name=category_name,
        short_label=category_name,
        dimensions="dimX=samplename dimY=viewpoint dimA=summary_method",
        sortOrder="samplename=+ viewpoint=+ summary_method=+",
        tracktype="bigWig",
        visibility="hide",
        dragAndDrop="subTracks",
        allButtonPair="off",
    )

    # Only add a group if this is an assembly hub
    if groups_file:
        composite.add_params(group=hub_name)

    composite.add_subgroups([subgroup_vp, subgroup_sample, subgroup_method])
    # composite.add_params(html=os.path.basename(stats_report))

    for bw in df.itertuples():
        t = trackhub.Track(
            name=f'{bw.samplename}_{bw.viewpoint}_{bw.method.replace("-summary", "")}',
            source=bw.fn,
            autoScale="off",
            tracktype="bigWig",
            windowingFunction="maximum",
            subgroups={
                "viewpoint": bw.viewpoint.lower(),
                "samplename": bw.samplename.lower(),
                "summary_method": bw.method.split("-")[0],
            },
            color=",".join([str(int(x * 255)) for x in color_mapping[bw.samplename]]),
        )

        # Only add a group if this is an assembly hub
        if groups_file:
            t.add_params(group=hub_name)

        composite.add_subtrack(t)

    trackdb.add_tracks(composite)

# Add viewpoints
t = trackhub.Track(
    name=os.path.basename(viewpoints).replace(".bigBed", "").capitalize(),
    source=viewpoints,
    tracktype="bigBed",
)

if genomes_file:
    t.add_params(group=hub_name)

trackdb.add_tracks(t)


##############
#  Stage hub #
##############

staging_tmp_dir = "hub_tmp_dir"

# Stage the hub
trackhub.upload.stage_hub(hub=hub, staging=staging_tmp_dir)

# Edit the hub.txt file to include the stats report as descriptionUrl
with open(os.path.join(staging_tmp_dir, f"{hub_name}.hub.txt"), "a") as hubtxt:
    hubtxt.write("\n")
    hubtxt.write(
        f'descriptionUrl {snakemake.config["genome"]["name"]}/{os.path.basename(capcruncher_report)}\n'
    )

# Copy to the new location
shutil.copytree(
    "hub_tmp_dir",
    snakemake.output[0],
    dirs_exist_ok=True,
)

# Delete the staged hub
shutil.rmtree("hub_tmp_dir")

# Copy the stats report to the correct location
shutil.copy(
    capcruncher_report,
    os.path.join(
        snakemake.config["hub"]["dir"],
        snakemake.config["genome"]["name"],
        os.path.basename(capcruncher_report),
    ),
)
