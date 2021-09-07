use human_panic::setup_panic;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyDict};
use pyo3::wrap_pyfunction;
use std::collections::HashMap;
use std::fs::File;
use std::io;

mod count_fragments;
mod fastq_deduplication;
mod utils;

#[pyfunction]
fn load_bincode(py: Python, path: String) -> PyResult<&PyDict> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let deserialised: HashMap<u64, u64> = bincode::deserialize_from(reader).unwrap();
    let deserialised_dict = &deserialised.into_py_dict(py);

    Ok(deserialised_dict)
}

#[pyfunction]
#[pyo3(name = "fastq_parse")]
#[pyo3(text_signature = "(fastq_files: List, output: str, /)")]
fn fastq_parse_py(fastq_files: Vec<String>, parsed_output: String) -> PyResult<String> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    fastq_deduplication::parse_fastqs(fastq_files, parsed_output.clone()).unwrap();
    Ok(parsed_output.to_string())
}

#[pyfunction]
#[pyo3(name = "fastq_find_duplicates")]
#[pyo3(text_signature = "(fastq_parsed_files: List, output: str, /)")]
fn fastq_find_duplicates_py(infiles: Vec<String>, outfile: String) -> PyResult<String> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    fastq_deduplication::identify_duplicates(&mut infiles.to_owned(), &outfile).unwrap();
    Ok(outfile)
}

#[pyfunction]
#[pyo3(name = "fastq_remove_duplicates")]
#[pyo3(text_signature = "(fastq_files: list, duplicates: str, outputs: str, /) -> dict")]
fn fastq_remove_duplicates_py(
    py: Python,
    fastq_files: Vec<String>,
    duplicates: String,
    outfiles: Vec<String>,
) -> PyResult<&PyDict> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    let stats = fastq_deduplication::remove_duplicates(fastq_files, duplicates, outfiles);
    Ok(stats.unwrap().into_py_dict(py))
}

//Function groups all slices by the parent id and counts the occurences of each restriction fragment combination.
#[pyfunction]
#[pyo3(name = "count_restriction_fragment_combinations")]
#[pyo3(text_signature = "(infile: str, outfile: str, remove_viewpoint: bool, n_threads: int, chunksize: int)")]
fn count_restriction_fragment_combinations_py(
    _py: Python,
    infile: String,
    outfile: String,
    remove_viewpoint: bool,
    n_threads: Option<usize>,
    chunksize: Option<usize>,
) -> PyResult<String> {
    
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    
    let counts = count_fragments::count_restriction_fragment_combinations(
        infile,
        chunksize,
        n_threads,
        remove_viewpoint,
    ).unwrap();

    count_fragments::restriction_fragment_counts_to_tsv(outfile.clone(), counts)?;

    Ok(outfile)
}

#[pymodule]
#[pyo3(name = "libcapcruncher")]
fn libcapcruncher(_py: Python, module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(fastq_parse_py, module)?)?;
    module.add_function(wrap_pyfunction!(fastq_find_duplicates_py, module)?)?;
    module.add_function(wrap_pyfunction!(fastq_remove_duplicates_py, module)?)?;
    module.add_function(wrap_pyfunction!(load_bincode, module)?)?;
    module.add_function(wrap_pyfunction!(
        count_restriction_fragment_combinations_py,
        module
    )?)?;

    Ok(())
}
