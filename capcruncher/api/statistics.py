from pydantic import BaseModel, computed_field
from typing import List, Optional, Union, Dict, TypeVar, Generic
import pathlib
import pandas as pd
from enum import Enum, IntEnum


class FastqDeduplicationStatistics(BaseModel):
    """Statistics for Fastq deduplication"""

    id: str = "unknown_sample"
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

    id: str = "unknown_sample"
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
            id=self.id,
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
        return pd.DataFrame(self.hist.items(), columns=[name, "count"]).assign(
            **{"read_number": read_number}
        ).sort_values(by=["count", name])

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
    

    def __add__(self, other: Union['ReadPairStat[int]', 'ReadPairStat[Histogram]', 'ReadPairStat[SliceNumberStats]']):
        
        read_1 = self.read1 + other.read1
        read_2 = self.read2 + other.read2 if self.read2 is not None else None
        
        instance_type = type(self.read1)
        
        rps =  ReadPairStat[instance_type](read1=read_1, read2=read_2)
        
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
