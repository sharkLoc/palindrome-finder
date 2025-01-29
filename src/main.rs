pub mod command_line;
pub mod exact_matches;
pub mod fasta_parsing;
pub mod output;
pub mod run_algorithm;
pub mod wfa;
mod adapters;

use anyhow::{Ok, Result};
use clap::Parser;
use command_line::PalinArgs;
use run_algorithm::run;
use std::time::Instant;

fn main() -> Result<()> {
    let global_timer = Instant::now();
    
    let args = PalinArgs::parse();

    run(&args)?;
    
    let elapsed = global_timer.elapsed();
    
    println!("Total elapsed time: {:.2?}", elapsed);
    println!();
    println!("---Settings---\n{}", args);

    Ok(())
}
