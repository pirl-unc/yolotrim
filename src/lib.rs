use pyo3::prelude::*;
use num_cpus;

use std::fmt; 
use flate2::write;
use flate2::Compression;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str;
use std::ops;
use regex::bytes::Regex;
use fastq::{parse_path, Record};
use std::cmp::*;




fn print_fastq_impl(input_filename : String) -> () {
    parse_path(Some(input_filename), |parser| {
        let ncpus : usize = num_cpus::get();
        let nthreads = if ncpus > 1 { ncpus - 1 } else { 1 };
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
            // we can initialize thread local variables here.
    
            for record_set in record_sets {
                for record in record_set.iter() {
                    match (
                            str::from_utf8(record.head()), 
                            str::from_utf8(record.seq()), 
                            str::from_utf8(record.qual())) {
                        (Ok(head), Ok(seq), Ok(qual)) => {
                            println!(
                                ">{name:?}\n{seq:?}\n+\n{qual:?}", 
                                name=head, 
                                seq=seq,
                                qual=qual
                            );
                        },
                        _ => {
                            println!("Error parsing record");
                        }
                    }
                    
                }
            }
            // The values we return (it can be any type implementing `Send`)
            // are collected from the different threads by
            // `parser.parallel_each` and returned. See doc for a description of
            // error handling.
            0
        }).expect("Invalid FASTQ file");
        println!("{}", results.iter().sum::<usize>());
    }).expect("Invalid compression");
}

#[derive(Debug)]
struct MatchStats { 
    num_total : u32,
    num_forwards_matches : u32, 
    num_backwards_matches : u32,
    num_double_matches : u32,
    num_too_short : u32,
}

#[derive(Debug)]
struct Entry {
    id : String,
    seq : String,
    qual : String,
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


pub fn writer(filename: &str) -> Box<dyn Write> {
    /*
    Copied from: https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561

    Write normal or compressed files seamlessly
    Uses the presence of a `.gz` extension to decide
    Attempting to have a file writer too
    */
    let path = Path::new(filename);
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    }
}


/// Implementation is taken from https://doi.org/10.1101/082214
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}


fn trim_fastq_impl(input_filename : String, output_filename: String, max_primer_length : usize, min_poly_a_length : usize, max_poly_a_length : usize) -> () {
    parse_path(Some(input_filename), |parser| {
        let ncpus : usize = num_cpus::get();
        let nthreads = if ncpus > 1 { ncpus - 1 } else { 1 };
        let max_trim_length = max_primer_length + max_poly_a_length;
        let min_read_length = 2 * max_primer_length + min_poly_a_length;

        let f_regex_str = format!("(A{{{}, {}}}([CGT]A{{{}, {}}})*)[ACGT]{{0, {}}}$", min_poly_a_length, max_poly_a_length, min_poly_a_length, max_poly_a_length, max_primer_length);
        let b_regex_str = format!("^[ACGT]{{0, {}}}((T{{{}, {}}}[ACG])*T{{{}, {}}})", max_primer_length,  min_poly_a_length, max_poly_a_length, min_poly_a_length, max_poly_a_length);

        let f_regex : Regex = Regex::new(&f_regex_str).unwrap();
        let b_regex : Regex = Regex::new(&b_regex_str).unwrap();

        let results: Vec<(MatchStats, Vec<Entry>)> = parser.parallel_each(nthreads, move |record_sets| {
            let mut too_short : u32 = 0;
            let mut total_records = 0;
            let mut forwards_matches = 0;
            let mut backwards_matches = 0;
            let mut double_matches = 0;
            let mut trimmed_records : Vec<Entry> = Vec::new();

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

                            trimmed_records.push(Entry{
                                id: str::from_utf8(id).unwrap().to_string(), 
                                seq: str::from_utf8(new_seq).unwrap().to_string(),
                                qual: str::from_utf8(new_quals).unwrap().to_string()});

                        },
                        (None, Some(b_match)) => {
                            backwards_matches += 1;

                            let n_trimmed = b_match.len();
                            
                            let new_seq = revcomp(&seq[n_trimmed..]);
                            let new_quals: Vec<u8> = (&qual[n_trimmed..]).iter().rev().map(|x| *x).collect();
                            
                            trimmed_records.push(Entry{
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
            combined_records.extend(entries);
        }
        println!("=== {total_match:?}/{total_records:?} (missed={missed:?}, double={double:?}, too_short={too_short:?})", 
            total_match=forwards_matches + backwards_matches, 
            total_records=total_records,
            missed=total_records - (forwards_matches + backwards_matches + double_matches),
            double=double_matches,
            too_short=too_short);
        let mut output_writer = writer(&output_filename);
        println!("Writing output to {output_filename}", output_filename=output_filename);
        for entry in combined_records {
            output_writer.write(format!("{}\n{}\n{}\n", entry.id, entry.seq, entry.qual)).unwrap();
    
        };
        output_writer.flush().unwrap();

    }).expect("Invalid compression");
}

#[pyfunction]
fn trim_fastq(input_filename : String, output_filename: String) -> PyResult<()> {
    trim_fastq_impl(input_filename, output_filename, 100, 10, 1000);
    Ok(())
}

#[pyfunction]
fn print_fastq(input_filename : String) -> PyResult<()> {
    print_fastq_impl(input_filename);
    Ok(())
}

#[pymodule]
fn yolotrim(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(trim_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(print_fastq, m)?)?;
    Ok(())
}