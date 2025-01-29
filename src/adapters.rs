use std::{fs::File, io::BufReader};

use anyhow::{Ok, Result};
use block_aligner::{
    cigar::Cigar, percent_len, scan_block::{AlignResult, Block, PaddedBytes}, scores::{Gaps, NucMatrix, NW1}
};

use crate::{command_line::AdapterArgs, fasta_parsing::{Fasta, FastaIterator}};

pub struct Adapter {
    pub name: String,
    pub ref_name: String,
    pub alignment: String,
    pub result: AlignResult,
}

impl Adapter {
    pub fn new(name: String, ref_name: String, alignment: String, result: AlignResult) -> Self {
        Self {
            name,
            ref_name,
            alignment,
            result
        }
    }
    pub fn get_name(&self) -> &str {
        &self.name
    }
    pub fn get_ref(&self) -> &str {
        &self.ref_name
    }
    pub fn get_seq(&self) -> &str {
        &self.alignment
    }
    pub fn get_result(&self) -> &AlignResult{
        &self.result
    }
}

pub fn align_adapters(sequence: Fasta, cmds: &AdapterArgs, output: &mut Vec<Adapter>) -> Result<()> {
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };

    let seq = sequence.get_sequence();
    let len = seq.len();

    let min_block_size: usize = percent_len(cmds.largest_adapter_len, 1.0);
    let max_block_size: usize = 8*min_block_size;

    let adapter_file = File::open(&cmds.adapters_file_path)?;

    let iterator = FastaIterator::new(BufReader::new(adapter_file));

    let read_padded = PaddedBytes::from_bytes::<NucMatrix>(seq.as_bytes(), max_block_size);

    let mut local = Block::<true, false, false, true, true>::new(min_block_size, len, max_block_size);

    for line in iterator {
        let line = line?;
        let adapter = line.get_sequence();

        let mut cigar = Cigar::new(min_block_size, len);
        
        
        let adapter_padded = PaddedBytes::from_bytes::<NucMatrix>(adapter.as_bytes(), max_block_size);
        local.align(
            &adapter_padded,
            &read_padded,
            &NW1,
            gaps,
            min_block_size..=max_block_size,
            0,
        );
        let res = local.res();
        local.trace().cigar_eq(&adapter_padded, &read_padded, res.query_idx, res.reference_idx, &mut cigar);
        let result = Adapter::new(
            sequence.get_name().split_once(' ').unwrap().0.to_owned(),
            line.get_name().to_owned(),
            cigar.to_string(),
            res
        );
        output.push(result);
    }
    Ok(())
}


