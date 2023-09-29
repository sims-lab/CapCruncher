from pydantic import BaseModel, computed_field
from typing import List, Optional, Union
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
    def percentage(self):
        return self.duplicates / self.total * 100
    
    @computed_field
    @property
    def unique(self):
        return self.total - self.duplicates

class FastqTrimmingStatistics(BaseModel):
    """Statistics for Fastq trimming"""
    id: str = "unknown_sample"
    reads_input: int
    reads_output: int
    reads_with_adapter_identified: int
    
    @computed_field
    @property
    def percentage_trimmed(self):
        return self.reads_with_adapter_identified / self.reads_input * 100

    @computed_field
    @property
    def percentage_passing_quality_filter(self):
        return self.reads_output / self.reads_input * 100
    
    @classmethod
    def from_multiqc_entry(cls, entry: pd.Series):
        cls(
            id=entry["Sample"],
            reads_input=entry["r_processed"],
            reads_output=entry["r_written"],
            reads_with_adapter_identified=entry["r_with_adapters"]
        )
    
    def __add__(self, other: "FastqTrimmingStatistics"):
        return FastqTrimmingStatistics(
            id=self.id,
            reads_input=self.reads_input + other.reads_input,
            reads_output=self.reads_output + other.reads_output,
            reads_with_adapter_identified=self.reads_with_adapter_identified + other.reads_with_adapter_identified
        )
    

class ReadNumber(IntEnum):
    """Enum for read number"""
    read1 = 1
    read2 = 2
    unpaired = 0
    

    
class HistogramLength(BaseModel):
    length: List[int]
    count: List[int]
    read_in_pair: ReadNumber
    
class HistogramNumber(BaseModel):
    statistic_name: str
    number: List[int]
    count: List[int]
    read_in_pair: ReadNumber
    
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, statistic_name: str, read_in_pair: ReadNumber):
        return cls(
            statistic_name=statistic_name,
            number=df[statistic_name].tolist(),
            count=df["count"].tolist(),
            read_in_pair=read_in_pair
        )
    
class FastqDigestionStatistics(BaseModel):
    """Statistics for Fastq digestion"""
    id: str = "unknown_sample"
    read_pairs_input: int
    read_pairs_output: int
    slice_lengths: HistogramLength
    slices_unfiltered: HistogramNumber
    slices_filtered: HistogramNumber
    
    @computed_field
    @property
    def percentage(self):
        return self.read_pairs_output / self.read_pairs_input * 100
    
    @computed_field
    @property
    def read_pairs_filtered(self):
        return self.read_pairs_input - self.read_pairs_output
    
    @computed_field
    @property
    def histogram_length(self):
        return pd.DataFrame([x.model_dump() for x in self.slice_lengths])

    @computed_field
    @property
    def histogram_unfiltered(self):
        return pd.DataFrame([x.model_dump() for x in self.slices_unfiltered])
    
    @computed_field
    @property
    def histogram_filtered(self):
        return pd.DataFrame([x.model_dump() for x in self.slices_filtered])
    
    
    
    
    


        
        
        
    
    
    
    
    
    

    
    