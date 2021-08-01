use pyo3::{prelude::*};
use pyo3::types::{PyDict, IntoPyDict};
use std::fs::File;
use std::io;
use std::collections::HashMap;
pub mod fastq_deduplication;


#[pyfunction]
fn load_bincode(py: Python, path: String) -> PyResult<&PyDict>{
    
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let deserialised: HashMap<u64, u64> = bincode::deserialize_from(reader).unwrap();
    let deserialised_dict = &deserialised.into_py_dict(py);

    Ok(deserialised_dict)

}

#[pyfunction]
#[pyo3(name = "fastq_parse")]
#[pyo3(text_signature = "(fastq_files:List, output: str, /)")]
fn fastq_parse_py(fastq_files: Vec<String>, parsed_output: String) -> PyResult<String> {
    
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    fastq_deduplication::parse_fastqs(fastq_files, parsed_output.clone()).unwrap();
    Ok(parsed_output.to_string())
}

#[pyfunction]
#[pyo3(name = "fastq_find_duplicates")]
#[pyo3(text_signature = "(fastq_parsed_files: List, output: str, /)")]
fn fastq_find_duplicates_py(json_input: Vec<String>, json_output: String) -> PyResult<String> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    fastq_deduplication::identify_duplicates(&mut json_input.to_owned(), &json_output).unwrap();
    Ok(json_output)
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

#[pymodule]
#[pyo3(name = "libcapcruncher")]
fn libcapcruncher(_py: Python, module: &PyModule) -> PyResult<()> {

    module.add_function(wrap_pyfunction!(fastq_parse_py, module)?)?;
    module.add_function(wrap_pyfunction!(fastq_find_duplicates_py, module)?)?;
    module.add_function(wrap_pyfunction!(fastq_remove_duplicates_py, module)?)?;
    module.add_function(wrap_pyfunction!(load_bincode, module)?)?;

    Ok(())
}
