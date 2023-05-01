use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str;
use std::ops;
use std::cmp::*;

use pyo3::prelude::*;

use num_cpus;
use regex::bytes::Regex;
use fastq::{parse_path, Record};
use rustc_hash::FxHashMap;
use indicatif::ProgressBar;
use gzp::{par::compress::{ParCompress, ParCompressBuilder}, deflate::Gzip};
use rayon::prelude::*;


/// Implementation is taken from https://doi.org/10.1101/082214
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}

pub fn num_threads() -> usize {
    num_cpus::get()
}


#[derive(Debug)]
struct MatchStats { 
    num_total : u32,
    num_forwards_matches : u32, 
    num_backwards_matches : u32,
    num_double_matches : u32,
    num_too_short : u32,
}


impl ops::Add<MatchStats> for MatchStats {
    type Output = MatchStats;

    fn add(self, other: MatchStats) -> MatchStats {
        MatchStats {
            num_total: self.num_total + other.num_total, 
            num_forwards_matches: self.num_forwards_matches + other.num_forwards_matches, 
            num_backwards_matches: self.num_backwards_matches + other.num_backwards_matches, 
            num_double_matches: self.num_double_matches + other.num_double_matches,
            num_too_short: self.num_too_short + other.num_too_short,
        }
    }
}



#[derive(Debug)]
struct FastqEntry {
    id : String,
    seq : String,
    qual : String,
}



fn count_prefixes(entries: &Vec<FastqEntry>, prefix_length: usize) -> FxHashMap<&str, u32> {
    let n = entries.len();
    entries
        .par_chunks(n / (num_threads()))
        .map(|entries| {
            let mut hash_map = FxHashMap::default();
            for entry in entries {
                *hash_map.entry(&entry.seq[..prefix_length]).or_insert(0) += 1;
            }; 
            hash_map

        })
    .reduce_with(|mut m1, m2| {
        for (k, v) in m2 {
            *m1.entry(k).or_default() += v;
        }
        m1
    })
    .unwrap()
}

fn infer_5p_primer(min_primer_length : usize, max_primer_length : usize, fastq_entries : &Vec<FastqEntry>) -> (String, usize, u32) {
    let mut best_prefix_length = min_primer_length;
    let mut best_prefix : String = "".to_string();
    let mut best_prefix_count : u32 = 0;

    for prefix_length in min_primer_length..max_primer_length {
        let prefix_counts: FxHashMap<&str, u32> = count_prefixes(&fastq_entries, prefix_length);
        let mut prefix_counts_vec : Vec<(&&str, &u32)> = prefix_counts.iter().collect();
        prefix_counts_vec.sort_by(|a, b| b.1.cmp(a.1));
        let (curr_best_prefix, curr_best_prefix_count) = prefix_counts_vec[0].clone();
        if (*curr_best_prefix_count as f64 / fastq_entries.len() as f64) < 0.7 {
            break;
        } else {
            best_prefix = curr_best_prefix.clone().to_string();
            best_prefix_length = prefix_length;
            best_prefix_count = *curr_best_prefix_count;
        }
    }
    (best_prefix, best_prefix_length, best_prefix_count)
}

fn trim_5p_primer_and_format(prefix : String, fastq_entries : &Vec<FastqEntry>) -> Vec<String> {
    let prefix_one_shorter = prefix[1..prefix.len()-1].to_string();
    let prefix_two_shorter = prefix[2..prefix.len()-2].to_string();
    let n = prefix.len(); 
    fastq_entries.par_iter().map(|entry| {
        let mut trimmed_seq = entry.seq.clone();
        let mut trimmed_qual = entry.qual.clone();
        if entry.seq.starts_with(&prefix) {
            trimmed_seq = trimmed_seq[n..].to_string();
            trimmed_qual = trimmed_qual[n..].to_string();
        } else if entry.seq.starts_with(&prefix_one_shorter) {
            trimmed_seq = trimmed_seq[(n - 1)..].to_string();
            trimmed_qual = trimmed_qual[(n - 1)..].to_string();
        } else if entry.seq.starts_with(&prefix_two_shorter) {
            trimmed_seq = trimmed_seq[(n - 2)..].to_string();
            trimmed_qual = trimmed_qual[(n - 2)..].to_string();
        };
        FastqEntry{id: entry.id.clone(), seq: trimmed_seq, qual: trimmed_qual}
    }).map(|entry| {
        format!("@{}\n{}\n+\n{}", entry.id, entry.seq, entry.qual)
    }).collect()
}

pub fn create_writer(filename: &str) -> Box<dyn Write> {
    /*
    Copied from: https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561

    Write normal or compressed files seamlessly
    Uses the presence of a `.gz` extension to decide
    Attempting to have a file writer too
    */
    let path = Path::new(filename);
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.to_string()),
        Ok(file) => file,
    };
    let default_capacity = 8 * 1024 * 1024;
    let buffered_file_writer = BufWriter::with_capacity(default_capacity, file);
    if path.extension() == Some(OsStr::new("gz")) {
        // let compression_level = Compression::default();
        let parallel_encoder: ParCompress<Gzip> = ParCompressBuilder::new().from_writer(buffered_file_writer);
        let buffered_encoder = BufWriter::with_capacity(default_capacity, parallel_encoder);
        Box::new(buffered_encoder)
    } else {
        Box::new(buffered_file_writer)
    }
}





fn trim_fastq_impl(
        input_filename : String, 
        output_filename: String, 
        min_primer_length : usize,
        max_primer_length : usize, 
        min_poly_a_length : usize, 
        max_poly_a_length : usize) -> () {
    println!("Reading FASTQ file {}...", input_filename);

    parse_path(Some(input_filename), |parser| {
        let max_trim_length = max_primer_length + max_poly_a_length;
        let min_read_length = 2 * max_primer_length + min_poly_a_length;

        let f_regex_str = format!("(A{{{}, {}}}([CGT]A{{{}, {}}})*)[ACGT]{{0, {}}}$", min_poly_a_length, max_poly_a_length, min_poly_a_length, max_poly_a_length, max_primer_length);
        let b_regex_str = format!("^[ACGT]{{0, {}}}((T{{{}, {}}}[ACG])*T{{{}, {}}})", max_primer_length,  min_poly_a_length, max_poly_a_length, min_poly_a_length, max_poly_a_length);

        let f_regex : Regex = Regex::new(&f_regex_str).unwrap();
        let b_regex : Regex = Regex::new(&b_regex_str).unwrap();

        let results: Vec<(MatchStats, Vec<FastqEntry>)> = parser.parallel_each(
                num_threads(), 
                move |record_sets| {
            let mut too_short : u32 = 0;
            let mut total_records = 0;
            let mut forwards_matches = 0;
            let mut backwards_matches = 0;
            let mut double_matches = 0;
            let mut trimmed_records : Vec<FastqEntry> = Vec::new();

            for record_set in record_sets {
                for record in record_set.iter() {
                    let id = record.head();
                    let seq = record.seq();
                    let qual = record.qual();
                    let n = seq.len(); 
                    total_records += 1;
                    if n < min_read_length {
                        too_short += 1; 
                        continue;
                    }
                    let slice_length : usize = min(n, max_trim_length);
    
                    let prefix = &seq[..slice_length];
                    
                    let suffix: &[u8] = &seq[n - slice_length ..];

                    let maybe_f_match = f_regex.find(suffix);
                    let maybe_b_match = b_regex.find(prefix);

                    match (maybe_f_match, maybe_b_match) {
                        (None, None) => {},
                        (Some(_), Some(_)) => {
                            double_matches += 1; 
                        },
                        (Some(f_match), None) => {
                            forwards_matches += 1;
                            
                            let n_trimmed = f_match.len();
                            let new_seq = &seq[..n - n_trimmed];
                            let new_quals = &qual[..n - n_trimmed];

                            trimmed_records.push(FastqEntry{
                                id: str::from_utf8(id).unwrap().to_string(), 
                                seq: str::from_utf8(new_seq).unwrap().to_string(),
                                qual: str::from_utf8(new_quals).unwrap().to_string()});

                        },
                        (None, Some(b_match)) => {
                            backwards_matches += 1;

                            let n_trimmed = b_match.len();
                            
                            let new_seq = revcomp(&seq[n_trimmed..]);
                            let new_quals: Vec<u8> = (&qual[n_trimmed..]).iter().rev().map(|x| *x).collect();
                            
                            trimmed_records.push(FastqEntry{
                                id: str::from_utf8(id).unwrap().to_string(), 
                                seq: str::from_utf8(&new_seq).unwrap().to_string(),
                                qual: str::from_utf8(&new_quals).unwrap().to_string()});
                        },
                    }                 
                }
                
            }
            let match_stats = MatchStats {
                num_forwards_matches: forwards_matches, 
                num_backwards_matches: backwards_matches,
                num_double_matches: double_matches,
                num_total: total_records,
                num_too_short: too_short,
            };
           (match_stats, trimmed_records)

        }).expect("Invalid FASTQ file");
        
        let mut total_records = 0;
        let mut forwards_matches = 0;
        let mut backwards_matches = 0;
        let mut double_matches = 0;
        let mut too_short : u32 = 0;
        let mut combined_records = Vec::new();

        for (match_stats, entries) in results {
            forwards_matches += match_stats.num_forwards_matches;
            backwards_matches += match_stats.num_backwards_matches; 
            double_matches += match_stats.num_double_matches;
            total_records += match_stats.num_total;
            too_short += match_stats.num_too_short;
            for entry in entries {
                if entry.seq.len() < max_primer_length {
                    too_short += 1;
                } else {
                    combined_records.push(entry);
                }
            }
        }
        println!("Finished 3' trimming: {total_match:?}/{total_records:?} reads (missing polyA = {missed:?}, 3'-3' artifacts = {double:?}, too short = {too_short:?})", 
            total_match=forwards_matches + backwards_matches, 
            total_records=total_records,
            missed=total_records - (forwards_matches + backwards_matches + double_matches),
            double=double_matches,
            too_short=too_short);
        
        println!("Inferring best 5' primer sequence...");
        let (best_prefix, best_prefix_length, best_prefix_count) = infer_5p_primer(
                min_primer_length, 
                max_primer_length, 
                &combined_records); 
        println!("Best 5' primer prefix: {} (length = {}, error free count = {})", best_prefix, best_prefix_length, best_prefix_count);
        println!("Trimming 5' ends of sequences and formatting...");
        let trimmed_strings = trim_5p_primer_and_format(best_prefix, &combined_records);
        println!("Writing output to {output_filename}", output_filename=output_filename);
        let mut output_writer = create_writer(&output_filename);
        let bar = ProgressBar::new(combined_records.len() as u64);
        for s in trimmed_strings {
            output_writer.write_all(s.as_bytes()).unwrap();
            bar.inc(1);
        }
        bar.finish();

    }).expect("Invalid compression");
}

#[pyfunction]
fn trim_fastq(
        input_filename : String,
        output_filename: String) -> PyResult<()> {
    trim_fastq_impl(input_filename, output_filename, 15, 70, 10, 1000);
    Ok(())
}

#[pymodule]
fn yolotrim(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(trim_fastq, m)?)?;
    Ok(())
}