from typing import Generic, TypeVar

import pandas as pd
from pydantic import BaseModel, computed_field, field_validator
from capcruncher.types import CisOrTrans, ReadType


def _normalise_read_type(value: str | ReadType) -> ReadType:
    if isinstance(value, ReadType):
        return value
    normalised = value.lower()
    if normalised == "unflashed":
        normalised = ReadType.PE
    return ReadType(normalised)


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
    read_number: int | float
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
    def from_multiqc_entry(cls, entry: pd.Series) -> "FastqTrimmingStatistics":
        return cls(
            sample=entry["sample"],
            read_number=entry["read_number"],
            reads_input=entry["r_processed"],
            reads_output=entry["r_written"],
            reads_with_adapter_identified=entry["r_with_adapters"],
        )

    def __add__(self, other: "FastqTrimmingStatistics") -> "FastqTrimmingStatistics":
        return FastqTrimmingStatistics(
            sample=self.sample,
            read_number=self.read_number,
            reads_input=self.reads_input + other.reads_input,
            reads_output=self.reads_output + other.reads_output,
            reads_with_adapter_identified=self.reads_with_adapter_identified
            + other.reads_with_adapter_identified,
        )


V = TypeVar("V")


class SliceNumberStats(BaseModel):
    unfiltered: int
    filtered: int

    def __add__(self, other: "SliceNumberStats") -> "SliceNumberStats":
        return SliceNumberStats(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
        )


class Histogram(BaseModel):
    name: str
    hist: dict[int, int]

    def to_dataframe(self, name: str = "value", read_number: str | None = None) -> pd.DataFrame:
        return (
            pd.DataFrame(self.hist.items(), columns=[name, "count"])
            .assign(**{"read_number": read_number})
            .sort_values(by=["count", name])
        )

    def __add__(self, other: "Histogram") -> "Histogram":
        return Histogram(
            name=self.name,
            hist={
                k: self.hist.get(k, 0) + other.hist.get(k, 0)
                for k in set(self.hist) | set(other.hist)
            },
        )


class ReadPairStat(BaseModel, Generic[V]):
    read1: Histogram | SliceNumberStats | int
    read2: Histogram | SliceNumberStats | int | None = None

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
        other: "ReadPairStat[int] | ReadPairStat[Histogram] | ReadPairStat[SliceNumberStats]",
    ) -> "ReadPairStat":
        read_1 = self.read1 + other.read1
        read_2 = self.read2 + other.read2 if self.read2 is not None else None

        return ReadPairStat(read1=read_1, read2=read_2)


class DigestionReadPairStats(BaseModel):
    unfiltered: ReadPairStat[int]
    filtered: ReadPairStat[int]

    def __add__(self, other: "DigestionReadPairStats") -> "DigestionReadPairStats":
        return DigestionReadPairStats(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
        )


class DigestionHistograms(BaseModel):
    unfiltered: ReadPairStat[Histogram]
    filtered: ReadPairStat[Histogram]
    lengths: ReadPairStat[Histogram]

    def __add__(self, other: "DigestionHistograms") -> "DigestionHistograms":
        return DigestionHistograms(
            unfiltered=self.unfiltered + other.unfiltered,
            filtered=self.filtered + other.filtered,
            lengths=self.lengths + other.lengths,
        )


class DigestionStats(BaseModel):
    sample: str
    read_type: ReadType
    read_stats: DigestionReadPairStats
    slice_stats: SliceNumberStats
    histograms: DigestionHistograms

    @field_validator("read_type", mode="before")
    @classmethod
    def validate_read_type(cls, value: str | ReadType) -> ReadType:
        return _normalise_read_type(value)

    def __add__(self, other: "DigestionStats") -> "DigestionStats":
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
    def percentage_combined(self) -> float:
        return self.n_combined / self.n_total * 100


class FlashOverallStats(BaseModel):
    samples: list[FlashStats]

    @classmethod
    def from_multiqc(cls, multiqc_data: str | pd.DataFrame) -> "FlashOverallStats":
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
    read_type: ReadType

    @field_validator("read_type", mode="before")
    @classmethod
    def validate_read_type(cls, value: str | ReadType) -> ReadType:
        return _normalise_read_type(value)

    @classmethod
    def from_slice_stats_dataframe(
        cls, df: pd.DataFrame, stage: str, sample: str, read_type: ReadType
    ) -> "SliceFilterStats":
        return cls(
            sample=sample,
            stage=stage,
            read_type=read_type,
            n_fragments=df["unique_fragments"],
            n_slices=df["unique_slices"],
        )


class SliceFilterStatsList(BaseModel):
    stats: list[SliceFilterStats]

    @classmethod
    def from_list(cls, stats: list[SliceFilterStats]) -> "SliceFilterStatsList":
        return cls(stats=stats)
    
    
class AlignmentDeduplicationStats(BaseModel):
    sample: str
    read_type: ReadType
    n_total_reads: int
    n_unique_reads: int
    n_total_slices: int
    n_unique_slices: int

    @field_validator("read_type", mode="before")
    @classmethod
    def validate_read_type(cls, value: str | ReadType) -> ReadType:
        return _normalise_read_type(value)
    
    @computed_field
    @property
    def percentage_unique_reads(self) -> float:
        return self.n_unique_reads / self.n_total_reads * 100
    
    @computed_field
    @property
    def n_duplicate_reads(self) -> int:
        return self.n_total_reads - self.n_unique_reads
    
    @computed_field
    @property
    def percentage_duplicate_slices(self) -> float:
        return self.n_duplicate_slices / self.n_total_slices * 100
    
    @computed_field
    @property
    def n_duplicate_slices(self) -> int:
        return self.n_total_slices - self.n_unique_slices
    


class CisOrTransStat(BaseModel):
    sample: str
    read_type: ReadType
    viewpoint: str
    cis_or_trans: CisOrTrans
    count: int

    @field_validator("read_type", mode="before")
    @classmethod
    def validate_read_type(cls, value: str | ReadType) -> ReadType:
        return _normalise_read_type(value)


class CisOrTransStats(BaseModel):
    stats: list[CisOrTransStat]

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "CisOrTransStats":
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
        return cls(stats=stats)
