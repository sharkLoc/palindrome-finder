use std::fs::File;

use anyhow::{anyhow, ensure, Result};

use crate::{
    command_line::{
        AlgorithmType::{FixedMismatch, Wfa}, FixedArgs, PalinArgs, WfaArgs
    },
    adapters::align,
    exact_matches::fixed_match,
    fasta_parsing::Fasta,
    output::PalindromeData,
    wfa::wfa_palins,
};

pub fn run(args: &PalinArgs, iterator: Box<dyn Iterator<Item = Result<Fasta>>>) -> Result<Vec<PalindromeData>> {
    let mut palins = Vec::new();
    match &args.command {
        Wfa(cmds) => run_wfa(args, cmds, iterator, &mut palins)?,
        FixedMismatch(cmds) => run_fixed_match(args, cmds, iterator, &mut palins)?,
    }
    Ok(palins)
}

pub fn write_output(args: &PalinArgs){
    
}

pub fn run_adapters(file_path: String, iterator: Box<dyn Iterator<Item = Result<Fasta>>>) -> Result<()> {
    let adapter_file = File::open(file_path);
    let adapter_file = match adapter_file {
        Result::Ok(adapter_file) => adapter_file,
        Err(err) => return Err(anyhow!("Invalid file format: {err}"))
    };
    let mut adapters: Vec<Fasta> = Vec::new();
    for line in iterator {
        let line = match line{
            Result::Ok(line) => line,
            Err(err) => return Err(anyhow!("Invalid file format: {err}"))
        };
        align(line, &adapter_file, &mut adapters)?
    }
    Ok(())
}

pub fn run_wfa(
    args: &PalinArgs,
    wfa_args: &WfaArgs,
    iterator: Box<dyn Iterator<Item = Result<Fasta>>>,
    output: &mut Vec<PalindromeData>,
) -> Result<()> {
    
    ensure!(wfa_args.match_bonus > 0.0, "Match bonus not positive");
    ensure!(
        wfa_args.mismatch_penalty > 0.0,
        "Mismatch penalty not positive"
    );
    ensure!(
        0.0 < wfa_args.mismatch_proportion && wfa_args.mismatch_proportion < 1.0,
        "Mismatch-length ratio not between 0 and 1"
    );
    ensure!(wfa_args.x_drop > 0.0, "X-drop not positive");

    let mut palins = Vec::new();
    for line in iterator {
        let line = line?;
        wfa_palins(line, output, args, wfa_args)?;
        output.append(&mut palins);
        palins.clear();
    }

    Ok(())
}

pub fn run_fixed_match(
    args: &PalinArgs,
    cmds: &FixedArgs,
    iterator: Box<dyn Iterator<Item = Result<Fasta>>>,
    output: &mut Vec<PalindromeData>,
) -> Result<()> {
    let mut palins = Vec::new();
    for line in iterator {
        let line = line?;
        fixed_match(line, output, args, cmds)?;
        output.append(&mut palins);
        palins.clear();
    }

    Ok(())
}
