
use anyhow::Result;

use crate::{
    adapters::align_adapters, command_line::{
        AdapterArgs, AlgorithmType::{FixedMismatch, Wfa, Adapters}, PalinArgs
    }, exact_matches::fixed_match, fasta_parsing::{parse, Fasta}, output::{write_adapters, write_palins, PalindromeData}, wfa::wfa_palins
};

pub fn run(args: &PalinArgs) -> Result<()> {
    let iterator = parse(args)?;
    let output_file = &args.output_file;
    
    match &args.command {
        Wfa(cmds) => run_algorithm(iterator, output_file, |fasta, palins| wfa_palins(fasta, palins, cmds))?,
        FixedMismatch(cmds) => run_algorithm(iterator, output_file, |fasta, palins| fixed_match(fasta, palins, cmds))?,
        Adapters(cmds) => run_adapters(cmds, iterator, output_file)?,
    }
    
    Ok(())
}

fn run_adapters(cmds: &AdapterArgs, iterator: Box<dyn Iterator<Item = Result<Fasta>>>, output_file: &str) -> Result<()> {
    let mut adapters = Vec::new();
    for fasta in iterator {
        align_adapters(fasta?, cmds, &mut adapters)?;
    }
    write_adapters(&mut adapters, output_file)?;
    Ok(())
}

fn run_algorithm<F>(iterator: Box<dyn Iterator<Item = Result<Fasta>>>, output_file: &str, algo: F) -> Result<()>
where
    F: Fn(Fasta, &mut Vec<PalindromeData>) -> Result<()>,
{
    let mut palins = Vec::new();
    for fasta in iterator {
        algo(fasta?, &mut palins)?;
    }

    write_palins(&mut palins, output_file)?;

    Ok(())
}
