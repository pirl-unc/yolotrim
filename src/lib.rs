use pyo3::prelude::*;
use num_cpus;
use fastq::{parse_path, Record};


fn trim_fastq_impl(input_filename : String, output_filename: String) -> () {
 
    parse_path(Some(input_filename), |parser| {
        let ncpus : usize = num_cpus::get();
        let nthreads = if ncpus > 1 { ncpus - 1 } else { 1 };
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
            // we can initialize thread local variables here.
    
            for record_set in record_sets {
                for record in record_set.iter() {
     
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

/// Formats the sum of two numbers as string.
#[pyfunction]
fn trim_fastq(input_filename : String, output_filename: String) -> PyResult<()> {
    trim_fastq_impl(input_filename, output_filename);
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn yolotrim(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(trim_fastq, m)?)?;
    Ok(())
}