use pyo3::prelude::*;
use num_cpus;
use std::str;
use regex::bytes::Regex;
use fastq::{parse_path, Record};
use lazy_static::lazy_static;


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

fn trim_fastq_impl(input_filename : String, output_filename: String) -> () {
 
    parse_path(Some(input_filename), |parser| {
        let ncpus : usize = num_cpus::get();
        let nthreads = if ncpus > 1 { ncpus - 1 } else { 1 };
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
            // we can initialize thread local variables here.
            lazy_static! {
                static ref RE: Regex = Regex::new(r"^[ACGTN]+$").unwrap();

            };
            for record_set in record_sets {
                for record in record_set.iter() {
                    let id = record.head();
                    let seq = record.seq();
                    let qual = record.qual();
                    RE.is_match(seq);
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

#[pyfunction]
fn trim_fastq(input_filename : String, output_filename: String) -> PyResult<()> {
    trim_fastq_impl(input_filename, output_filename);
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