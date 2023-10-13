from pydantic import BaseModel, computed_field
from typing import List, Optional, Union, Dict, TypeVar, Generic
import pathlib
import pandas as pd
from enum import Enum, IntEnum
from typing import Literal


class ReadType(Enum):
    flashed: str = "flashed"
    pe: str = "unflashed"


class FastqDeduplicationStatistics(BaseModel):
    """Statistics for Fastq deduplication"""

    sample: str = "unknown_sample"
    total: int
    duplicates: int

    @computed_field
    @property
    def percentage(self) -> float:
        return self.duplicates / self.total * 100

    @computed_field
    @property
    def unique(self) -> int:
        return self.total - self.duplicates


class FastqTrimmingStatistics(BaseModel):
    """Statistics for Fastq trimming"""

    sample: str = "unknown_sample"
    reads_input: int
    reads_output: int
    reads_with_adapter_identified: int

    @computed_field
    @property
    def percentage_trimmed(self) -> float:
        return self.reads_with_adapter_identified / self.reads_input * 100

    @computed_field
    @property
    def percentage_passing_quality_filter(self) -> float:
        return self.reads_output / self.reads_input * 100

    @classmethod
    def from_multiqc_entry(cls, entry: pd.Series):
        cls(
            id=entry["Sample"],
            reads_input=entry["r_processed"],
            reads_output=entry["r_written"],
            reads_with_adapter_identified=entry["r_with_adapters"],
        )

    def __add__(self, other: "FastqTrimmingStatistics"):
        return FastqTrimmingStatistics(
            sample=self.sample,
            reads_input=self.reads_input + other.reads_input,
            reads_output=self.reads_output + other.reads_output,
            reads_with_adapter_identified=self.reads_with_adapter_identified
            + other.reads_with_adapter_identified,
        )


V = TypeVar('V')


class SliceNumberStats(BaseModel):
    unfiltered: int
    filtered: int

    def __add__(self, other: "SliceNumberStats"):
        return SliceNumberStats(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
        )


class Histogram(BaseModel):
    name: str
    hist: Dict[int, int]

    def to_dataframe(
        self, name: Optional[str] = "value", read_number: Optional[str] = None
    ):
        return (
            pd.DataFrame(self.hist.items(), columns=[name, "count"])
            .assign(**{"read_number": read_number})
            .sort_values(by=["count", name])
        )

    def __add__(self, other: "Histogram"):
        return Histogram(
            name=self.name,
            hist={
                k: self.hist.get(k, 0) + other.hist.get(k, 0)
                for k in set(self.hist) | set(other.hist)
            },
        )


class ReadPairStat(BaseModel, Generic[V]):
    read1: Union[Histogram, SliceNumberStats, int]
    read2: Optional[Union[Histogram, SliceNumberStats, int]] = None

    def to_dataframe(self) -> pd.DataFrame:
        frames = []
        frames.append(
            self.read1.to_dataframe(read_number="read1", name=self.read1.name)
        )
        if self.read2 is not None:
            frames.append(
                self.read2.to_dataframe(read_number="read2", name=self.read2.name)
            )

        return pd.concat(frames)

    def __add__(
        self,
        other: Union[
            'ReadPairStat[int]',
            'ReadPairStat[Histogram]',
            'ReadPairStat[SliceNumberStats]',
        ],
    ):
        read_1 = self.read1 + other.read1
        read_2 = self.read2 + other.read2 if self.read2 is not None else None

        instance_type = type(self.read1)

        rps = ReadPairStat[instance_type](read1=read_1, read2=read_2)

        return rps


class DigestionReadPairStats(BaseModel):
    unfiltered: ReadPairStat[int]
    filtered: ReadPairStat[int]

    def __add__(self, other: "DigestionReadPairStats"):
        return DigestionReadPairStats(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
        )


class DigestionHistograms(BaseModel):
    unfiltered: ReadPairStat[Histogram]
    filtered: ReadPairStat[Histogram]
    lengths: ReadPairStat[Histogram]

    def __add__(self, other):
        return DigestionHistograms(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
            lengths=self.lengths + other.lengths,
        )


class DigestionStats(BaseModel):
    sample: str
    read_type: str
    read_stats: DigestionReadPairStats
    slice_stats: SliceNumberStats
    histograms: DigestionHistograms

    def __add__(self, other) -> 'DigestionStats':
        return DigestionStats(
            sample=self.sample,
            read_type=self.read_type,
            histograms=self.histograms + other.histograms,
            read_stats=self.read_stats + other.read_stats,
            slice_stats=self.slice_stats + other.slice_stats,
        )


class FlashStats(BaseModel):
    sample: str
    n_combined: int
    n_uncombined: int

    @computed_field
    @property
    def n_total(self) -> int:
        return self.n_combined + self.n_uncombined

    @computed_field
    @property
    def percentage_combined(self) -> int:
        return self.n_combined / self.n_total * 100


class FlashOverallStats(BaseModel):
    samples: List[FlashStats]

    @classmethod
    def from_multiqc(cls, multiqc_data: Union[str, pd.DataFrame]):
        if isinstance(multiqc_data, str):
            multiqc_data = pd.read_csv(multiqc_data, sep="\t")

        multiqc_data["sample_name"] = multiqc_data["Sample"].str.split("_part").str[0]
        multiqc_data = (
            multiqc_data[["sample_name", "combopairs", "uncombopairs"]]
            .groupby("sample_name")
            .sum()
            .reset_index()
        )

        samples = [
            FlashStats(
                sample=row["sample_name"],
                n_combined=row["combopairs"],
                n_uncombined=row["uncombopairs"],
            )
            for _, row in multiqc_data.iterrows()
        ]

        return cls(samples=samples)


class SliceFilterStats(BaseModel):
    sample: str
    stage: str
    n_fragments: int
    n_slices: int
    read_type: str

    @classmethod
    def from_slice_stats_dataframe(
        cls, df: pd.DataFrame, stage: str, sample: str, read_type: str
    ):
        return cls(
            sample=sample,
            stage=stage,
            read_type=read_type,
            n_fragments=df["unique_fragments"],
            n_slices=df["unique_slices"],
        )


class SliceFilterStatsList(BaseModel):
    stats: List[SliceFilterStats]

    @classmethod
    def from_list(cls, stats: List[SliceFilterStats]):
        return cls(stats=stats)


class CisOrTransStat(BaseModel):
    sample: str
    read_type: str
    viewpoint: str
    cis_or_trans: Literal["cis", "trans"]
    count: int


class CisOrTransStats(BaseModel):
    stats: List[CisOrTransStat]

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame):
        stats = []
        for row in df.itertuples():
            stats.append(
                CisOrTransStat(
                    sample=row.sample,
                    read_type=row.read_type,
                    viewpoint=row.viewpoint,
                    cis_or_trans=row.cis_or_trans,
                    count=row.count,
                )
            )
