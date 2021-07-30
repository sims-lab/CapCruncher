use bincode;
use bio::io::fastq;
use flate2::{bufread, write, Compression};
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::path::Path;
use twox_hash::xxh3::hash64_with_seed;

fn get_reader_handle(path: &str) -> Box<dyn io::Read> {
    if path.ends_with(".gz") {
        let f = File::open(path).unwrap();
        Box::new(bufread::GzDecoder::new(io::BufReader::new(f)))
    } else {
        Box::new(File::open(path).unwrap())
    }
}

fn get_writer_handle(path: &str) -> Box<dyn io::Write> {
    let f = File::create(path).expect("Cannot open output file");
    if path.ends_with(".gz") {
        Box::new(io::BufWriter::new(write::GzEncoder::new(
            f,
            Compression::default(),
        )))
    } else {
        Box::new(io::BufWriter::new(f))
    }
}

struct FastqReaders<R: fastq::FastqRead> {
    readers: Vec<R>,
}

impl<R: fastq::FastqRead> FastqReaders<R> {
    pub fn new(readers: Vec<R>) -> Self {
        FastqReaders { readers }
    }
}

impl FastqReaders<fastq::Reader<Box<dyn io::Read>>> {
    pub fn from_files<P: AsRef<Path> + Debug>(paths: Vec<P>) -> Self {
        let mut readers = Vec::new();
        for p in paths {
            let path_str = p.as_ref().to_str().expect("Cannot convert to string");
            let file_handle = get_reader_handle(path_str);
            let reader = fastq::Reader::new(file_handle);
            readers.push(reader);
        }

        FastqReaders::new(readers)
    }
}

impl<R: fastq::FastqRead> Iterator for FastqReaders<R> {
    type Item = Vec<fastq::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let n_readers = self.readers.len();
        let mut records = Vec::with_capacity(n_readers);
        let mut all_records_ok = true;

        for i in 0..n_readers {
            let mut record = fastq::Record::new();
            self.readers[i]
                .read(&mut record)
                .expect("Error reading record");

            if !record.is_empty() {
                records.push(record);
            } else {
                all_records_ok = false;
            }
        }

        if all_records_ok {
            Some(records)
        } else {
            None
        }
    }
}

fn hash_records(records: Vec<fastq::Record>, number_of_records: usize) -> (u64, u64) {
    let mut ids = Vec::with_capacity(number_of_records);
    let mut seqs = Vec::with_capacity(number_of_records);

    for r in records {
        ids.push(r.id().to_owned());
        seqs.push(r.seq().to_owned());
    }

    let ids_hashed = hash64_with_seed(&ids.concat().as_bytes(), 42);
    let seqs_hashed = hash64_with_seed(&seqs.concat(), 42);
    (ids_hashed, seqs_hashed)

    // (ids.concat().to_string(), String::from_utf8(seqs.concat()).unwrap())
}

pub fn parse_fastqs<P: AsRef<Path> + Debug>(
    files: Vec<P>,
    output: P,
) -> Result<(), Box<dyn Error>> {
    let number_of_files = files.len();
    let fastqs = FastqReaders::from_files(files);
    let outfile = output.as_ref().to_owned();
    let mut hashed_details = HashMap::new();
    let writer = io::BufWriter::new(File::create(outfile)?);

    for (ii, records) in fastqs.enumerate() {
        let (id, seq) = hash_records(records, number_of_files);
        hashed_details.insert(id, seq);
        if ii % 10000 == 0 {
            println!("Processed {} reads", ii);
        }
    }

    bincode::serialize_into(writer, &hashed_details)?;
    Ok(())
}

pub fn identify_duplicates<P: AsRef<Path>>(
    parsed: &mut Vec<P>,
    outfile: &P,
) -> Result<(), Box<dyn Error>> {
    parsed.shuffle(&mut thread_rng());
    let mut unique_sequences: HashSet<u64> = HashSet::new();
    let mut duplicated_ids: HashSet<u64> = HashSet::new();

    for encoded in parsed {
        let file = File::open(encoded)?;
        let reader = std::io::BufReader::new(file);

        let parsed_fastq_data: HashMap<u64, u64> = bincode::deserialize_from(reader)?;

        for (hashed_id, hashed_seq) in parsed_fastq_data {
            if !unique_sequences.contains(&hashed_seq) {
                unique_sequences.insert(hashed_seq);
            } else {
                duplicated_ids.insert(hashed_id);
            }
        }
    }

    println!("Number of duplicates {}", duplicated_ids.len());

    let writer = io::BufWriter::new(File::create(outfile)?);
    bincode::serialize_into(writer, &duplicated_ids)?;

    Ok(())
}

pub fn remove_duplicates<P: AsRef<Path> + Debug>(
    infiles: Vec<P>,
    duplicates: P,
    outfiles: Vec<P>,
) -> Result<HashMap<String, u64>, Box<dyn Error>> {
    
    let number_of_files = infiles.len();
    let fastqs = FastqReaders::from_files(infiles);

    let duplicate_file = File::open(duplicates)?;
    let duplicate_reader = std::io::BufReader::new(duplicate_file);
    let duplicate_ids: HashSet<u64> = bincode::deserialize_from(duplicate_reader)?;

    let mut outfiles_handles: Vec<_> = outfiles
        .iter()
        .map(|f| fastq::Writer::new(get_writer_handle(f.as_ref().to_str().unwrap())))
        .collect();

    let mut counter_unique = 0;
    let mut counter_removed = 0;

    for (ii, records) in fastqs.enumerate() {
        if ii % 10000 == 0 {
            println!("{}", ii);
        }

        let id = records
            .iter()
            .map(|r| r.id().to_owned())
            .collect::<Vec<String>>();
        let id_hashed = hash64_with_seed(id.concat().as_bytes(), 42);
        
        if duplicate_ids.contains(&id_hashed) {
            counter_removed += 1;
        } else {
            counter_unique += 1;

            for ii in 0..number_of_files {
                outfiles_handles[ii].write_record(&records[ii])?;
            }
        }
    }

    let mut stats = HashMap::with_capacity(3);
    stats.insert("total".to_owned(), counter_unique + counter_removed);
    stats.insert("unique".to_owned(), counter_unique);
    stats.insert("removed".to_owned(), counter_removed);

    Ok(stats)
}

#[cfg(test)]
mod tests_mp {
    use super::*;

    // #[test]
    // fn test_parsing_mp() {
    //     // let fq1 =
    //     //     "/stopgap/milne_group/asmith/ccanalyser/capcruncher/data/test/duplicated_1.fastq.gz";
    //     // let fq2 =
    //     //     "/stopgap/milne_group/asmith/ccanalyser/capcruncher/data/test/duplicated_2.fastq.gz";

    //     let fq1 = "../test_rust_cc/Erythroid-1_1.fastq.gz";
    //     let fq2 = "../test_rust_cc/Erythroid-1_2.fastq.gz";

    //     let fq_files = vec![&fq1, &fq2];

    //     parse_fastqs(fq_files, "parsed_XX.json");
    // }
}
