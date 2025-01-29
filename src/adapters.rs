use std::{fs::File, io::BufReader};

use anyhow::{anyhow, Ok, Result};
use block_aligner::{scan_block::*, scores::*};

use crate::fasta_parsing::{Fasta, FastaIterator};


pub fn align(sequence: Fasta, adapter_file: &File, output: &mut Vec<Fasta>) -> Result<()> {
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };
    let min_block_size: usize = 32;
    let max_block_size: usize = 256;

    
    let iterator = FastaIterator::new(BufReader::new(adapter_file));

    let read_padded = PaddedBytes::from_bytes::<NucMatrix>(sequence.get_sequence().as_bytes(), max_block_size);

    for line in iterator {
        let line = match line {
            Result::Ok(line) => line,
            Err(err) => return Err(anyhow!("Invalid file format: {err}")),
        };
        let adapter = line.get_sequence();
        let adapter_padded =
            PaddedBytes::from_bytes::<NucMatrix>(adapter.as_bytes(), max_block_size);
        let a = Block::<_, true, false>::align(
            &adapter_padded,
            &read_padded,
            &NW1,
            gaps,
            min_block_size..=max_block_size,
            0,
        );
        let res = a.res();
        let result = Fasta::new(sequence.get_name().to_owned(), a.trace().cigar(res.query_idx, res.reference_idx).to_string());

        output.push(result);
    }
    Ok(())
}
