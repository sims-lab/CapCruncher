use human_panic::setup_panic;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyDict};
use pyo3::{wrap_pyfunction, wrap_pymodule};
use std::collections::HashMap;
use std::fs::File;
use std::io;

mod count_fragments;
mod fastq_deduplication;
mod utils;

use fastq_deduplication::*;
use count_fragments::*;


#[pymodule]
#[pyo3(name = "libcapcruncher")]
fn libcapcruncher(py: Python, module: &PyModule) -> PyResult<()> {
    module.add_wrapped(wrap_pymodule!(fastq_deduplication))?;
    module.add_wrapped(wrap_pymodule!(count_fragments))?;

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;
    sys_modules.set_item("libcapcruncher.fastq_deduplication", module.getattr("fastq_deduplication")?)?;
    sys_modules.set_item("libcapcruncher.count_fragments", module.getattr("count_fragments")?)?;

    Ok(())
}