use crate::utils;
use csv::ByteRecord;
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

//impl Iterator<Item = &'a DigestedReadRestrictionFragments<'a>>

fn count_restriction_fragment_combinations_in_chunk(
    slices: &mut Vec<DigestedReadRestrictionFragments>,
) -> HashMap<(i64, i64), i64> {
    let mut restriction_fragment_counts = HashMap::new();
    //slices.sort_unstable_by(|a, b| a.parent_read.cmp(&b.parent_read));

    let slices_groups = slices.iter().group_by(|&r| &r.parent_read);

    for (_parent_read, slice_group) in slices_groups.into_iter() {
        let restriction_fragments = slice_group
            .into_iter()
            .map(|s| s.restriction_fragment)
            .filter(|frag| *frag >= 0)
            .collect_vec();

        for comb in restriction_fragments.into_iter().combinations(2) {
            let mut f1 = comb[0];
            let mut f2 = comb[1];

            if f1 > f2 {
                std::mem::swap(&mut f1, &mut f2)
            }

            *restriction_fragment_counts
                .entry((f1, f2))
                .or_insert(0 as i64) += 1 as i64;
        }
    }

    restriction_fragment_counts
}

pub fn count_restriction_fragment_combinations<P: AsRef<Path>>(
    path: P,
    chunksize: Option<usize>,
    n_threads: Option<usize>,
    remove_viewpoint: bool,
) -> Result<BTreeMap<(i64, i64), i64>, Box<dyn Error>> {
    // Open TSV and read headers
    let mut reader = utils::get_tsv_reader(path)?;
    let headers = utils::get_tsv_headers(&mut reader)?;

    // Using a threadpool to process each chunk of the TSV file
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads.unwrap_or(4))
        .build()
        .unwrap();
    let (tx, rx) = channel();

    let tsv_chunks = reader
        .into_byte_records()
        .chunks(chunksize.unwrap_or(2e6 as usize));

    for chunk in &tsv_chunks {
        let slices_chunk = chunk.collect_vec();
        let mut slices_deserialised: Vec<DigestedReadRestrictionFragments> = slices_chunk
            .par_iter()
            .map(|res| res.as_ref().unwrap().deserialize(Some(&headers)).unwrap())
            .filter(|s: &DigestedReadRestrictionFragments| {
                s.capture == "." || remove_viewpoint == false
            })
            .collect();

        slices_deserialised.par_sort_by_key(|drf| drf.parent_read.clone());

        //println!("{:?}", &slices_deserialised);

        let tx = tx.clone();

        pool.spawn(move || {
            let counts = count_restriction_fragment_combinations_in_chunk(&mut slices_deserialised);
            tx.send(counts).expect("Cant send data")
        })
    }

    drop(tx);

    let mut restriction_fragment_counts_combined = BTreeMap::new();
    for counts in rx.iter() {
        //println!("{:?}", counts);
        for (k, v) in counts {
            *restriction_fragment_counts_combined.entry(k).or_insert(0) += v;
        }
    }

    Ok(restriction_fragment_counts_combined)
}

pub fn restriction_fragment_counts_to_tsv<P: AsRef<Path>>(
    path: P,
    counts: BTreeMap<(i64, i64), i64>,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counting() {
        let reporters =
            "/home/asmith/Projects/CapCruncher/capcruncher/data/test/Slc25A37_reporters.tsv.gz";
        let res = count_restriction_fragment_combinations(reporters, None, None, false);
        assert!(res.is_ok(), "Counting failed")
    }
}
