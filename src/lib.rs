use pyo3::prelude::*;
use num_cpus;
use std::str;
use std::ops;
use regex::bytes::Regex;
use fastq::{parse_path, Record};
use std::process;




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

fn trim_fastq_impl(input_filename : String, output_filename: String, max_primer_length : usize, min_poly_a_length : usize, max_poly_a_length : usize) -> () {
 
    parse_path(Some(input_filename), |parser| {
        let ncpus : usize = num_cpus::get();
        let nthreads = if ncpus > 1 { ncpus - 1 } else { 1 };
        let max_trim_length = max_primer_length + max_poly_a_length;
        let min_read_length = 2 * max_primer_length + max_poly_a_length;

        let b_regex_str = format!("^[ACGT]{{0, {}}}T{{{}, {}}}", max_primer_length, min_poly_a_length, max_poly_a_length);
        
        let f_regex_str = format!("A{{{}, {}}}[ACGT]{{0, {}}}$", min_poly_a_length, max_poly_a_length, max_primer_length);
        println!("{}", b_regex_str);
        println!("{}", f_regex_str);
        let b_regex : Regex = Regex::new(&b_regex_str).unwrap();
        let f_regex : Regex = Regex::new(&f_regex_str).unwrap();

        let results: Vec<MatchStats> = parser.parallel_each(nthreads, move |record_sets| {
            let mut too_short : u32 = 0;
            let mut total_records = 0;
            let mut forwards_matches = 0;
            let mut backwards_matches = 0;
            let mut double_matches = 0;
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
                    let forwards_match=f_regex.is_match(&seq[..max_trim_length]);
                    let backwards_match=b_regex.is_match(&seq[(n - max_trim_length) .. ]);
                    if forwards_match { 
                        if backwards_match {
                            double_matches += 1;
                        } else {
                            forwards_matches += 1; 
                        }
                    }
                    else if backwards_match { 
                        backwards_matches += 1; 
                    }
                 
                    /* 
                    
                    println!("f={forwards_match:?} b={backwards_match:?} | {seq_start:?}...{seq_end:?}", 
                        forwards_match=FORWARD_RE.is_match(seq),
                        backwards_match=BACKWARD_RE.is_match(seq),
                        seq_start=str::from_utf8(&seq[..100]).unwrap(), 
                        seq_end=str::from_utf8(&seq[n - 100..]).unwrap(), 
                    
                    );
                    */
                    
                }
                
                println!("{total_match:?}/{total_records:?} (missed={missed:?}, double={double:?}, too_short={too_short:?})", 
                    total_match=forwards_matches + backwards_matches, 
                    total_records=total_records,
                    missed=total_records - (forwards_matches + backwards_matches + double_matches),
                    double=double_matches,
                    too_short=too_short);
                
            }
            MatchStats {
                num_forwards_matches: forwards_matches, 
                num_backwards_matches: backwards_matches,
                num_double_matches: double_matches,
                num_total: total_records,
                num_too_short: too_short,
            }
        }).expect("Invalid FASTQ file");
        let mut total_records = 0;
        let mut forwards_matches = 0;
        let mut backwards_matches = 0;
        let mut double_matches = 0;
        let mut too_short : u32 = 0; 
        for  match_stats in results {
            forwards_matches += match_stats.num_forwards_matches;
            backwards_matches += match_stats.num_backwards_matches; 
            double_matches += match_stats.num_double_matches;
            total_records += match_stats.num_total;
            too_short += match_stats.num_too_short;
        }
        println!("=== {total_match:?}/{total_records:?} (missed={missed:?}, double={double:?}, too_short={too_short:?})", 
            total_match=forwards_matches + backwards_matches, 
            total_records=total_records,
            missed=total_records - (forwards_matches + backwards_matches + double_matches),
            double=double_matches,
            too_short=too_short);
    }).expect("Invalid compression");
}

#[pyfunction]
fn trim_fastq(input_filename : String, output_filename: String) -> PyResult<()> {
    trim_fastq_impl(input_filename, output_filename, 100, 10, 600);
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