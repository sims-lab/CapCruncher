use crate::utils;
use itertools::Itertools;
use pyo3::prelude::*;
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fmt::Debug;
use std::io;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc::channel;

#[derive(Debug, Deserialize, Clone)]
struct DigestedReadRestrictionFragments {
    parent_read: String,
    restriction_fragment: i64,
    capture: String,
}

fn count_restriction_fragment_combinations_in_chunk(
    slices: &mut Vec<DigestedReadRestrictionFragments>,
) -> HashMap<(i64, i64), usize> {
    let mut restriction_fragment_counts = HashMap::new();
    slices.sort_by_key(|r| r.parent_read.clone());
    let slices_groups = slices.iter().group_by(|r| &r.parent_read);

    for (_parent_read, slice_group) in slices_groups.into_iter() {
        let restriction_fragments = slice_group
            .into_iter()
            .map(|s| s.restriction_fragment)
            .filter(|frag| *frag >= 0);

        for comb in restriction_fragments.combinations(2) {
            let mut f1 = comb[0];
            let mut f2 = comb[1];

            if f1 > f2 {
                std::mem::swap(&mut f1, &mut f2)
            }

            restriction_fragment_counts
                .entry((f1, f2))
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }
    restriction_fragment_counts
}

fn get_tsv_reader<P: AsRef<Path>>(path: P) -> Result<csv::Reader<Box<dyn io::Read>>, io::Error> {
    let fh = utils::get_reader_handle(path.as_ref().to_str().expect("Not UTF-8 file name"))?;
    let reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(fh);
    Ok(reader)
}

fn get_tsv_headers<R: io::Read>(
    reader: &mut csv::Reader<R>,
) -> Result<csv::ByteRecord, Box<dyn Error>> {
    let headers = reader.headers().unwrap().as_byte_record().to_owned();
    Ok(headers)
}

pub fn count_restriction_fragment_combinations<P: AsRef<Path>>(
    path: P,
    chunksize: Option<usize>,
    n_threads: Option<usize>,
    remove_viewpoint: bool,
) -> Result<BTreeMap<(i64, i64), usize>, Box<dyn Error>> {
    // Open TSV and read headers
    let mut reader = get_tsv_reader(path)?;
    let headers = get_tsv_headers(&mut reader)?;

    // Using a threadpool to process each chunk of the TSV file
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads.unwrap_or(4))
        .build()
        .unwrap();
    let (tx, rx) = channel();

    // Run the counting operation on chunks of the TSV
    for chunk in &reader
        .byte_records()
        .chunks(chunksize.unwrap_or(2e6 as usize))
    {
        let mut slices = match remove_viewpoint {
            true => chunk
                .map(|r| {
                    r.unwrap()
                        .deserialize(Some(&headers))
                        .expect("Record does not match the expected structure")
                })
                .filter(|s: &DigestedReadRestrictionFragments| s.capture != ".")
                .collect(),
            false => chunk
                .map(|r| {
                    r.unwrap()
                        .deserialize(Some(&headers))
                        .expect("Record does not match the expected structure")
                })
                .collect(),
        };
        let tx = tx.clone();
        pool.spawn(move || {
            let counts = count_restriction_fragment_combinations_in_chunk(&mut slices);
            tx.send(counts).expect("Cant send data")
        })
    }

    drop(tx);

    // Aggregate the results
    let mut restriction_fragment_counts_combined = BTreeMap::new();
    for counts in rx.iter() {
        for (k, v) in counts {
            restriction_fragment_counts_combined
                .entry(k)
                .and_modify(|old| *old += v)
                .or_insert(v);
        }
    }

    Ok(restriction_fragment_counts_combined)
}

pub fn restriction_fragment_counts_to_tsv<P: AsRef<Path>>(
    path: P,
    counts: BTreeMap<(i64, i64), usize>,
) -> Result<(), io::Error> {
    let fh = utils::get_writer_handle(path.as_ref().to_str().expect("Not UTF-8 file name"))?;
    let mut writer = io::BufWriter::new(fh);

    writer
        .write_all(b"bin1_id\tbin2_id\tcount\n")
        .expect("Failed to write");
    for ((bin1, bin2), count) in counts {
        writer
            .write_all(format!("{}\t{}\t{}\n", bin1, bin2, count).as_bytes())
            .expect("Failed to write");
    }

    Ok(())
}

// Python bindings

/// Groups all slices by the parent id and counts the occurences of each restriction fragment combination.
#[pyfunction]
#[pyo3(name = "count_restriction_fragment_combinations")]
#[pyo3(
    text_signature = "(infile: str, outfile: str, remove_viewpoint: bool, n_threads: int, chunksize: int)"
)]
fn count_restriction_fragment_combinations_py(
    _py: Python,
    infile: String,
    outfile: String,
    remove_viewpoint: bool,
    n_threads: Option<usize>,
    chunksize: Option<usize>,
) -> PyResult<String> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    let counts =
        count_restriction_fragment_combinations(infile, chunksize, n_threads, remove_viewpoint)
            .unwrap();

    restriction_fragment_counts_to_tsv(outfile.clone(), counts)?;

    Ok(outfile)
}

#[pymodule]
#[pyo3(name = "count_fragments")]
fn fastq_deduplication(_py: Python, module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(
        count_restriction_fragment_combinations_py,
        module
    )?)?;
    Ok(())
}
