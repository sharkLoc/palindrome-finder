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

pub fn align_adapters(sequence: &mut Fasta, cmds: &AdapterArgs, output: &mut Vec<Adapter>) -> Result<()> {
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };

    let mut offset = 0;

    if cmds.remove_t{
        let t_res = remove_t(sequence);
        let offset = match t_res {
            Some(t_res) => t_res.reference_idx,
            None => 0
        };
    }
    let seq = sequence.get_sequence();
    let len = seq.len();

    let min_block_size: usize = percent_len(cmds.largest_adapter_len, 1.0);
    let max_block_size: usize = 8*min_block_size;

    let adapter_file = File::open(&cmds.adapters_file_path)?;

    let iterator = FastaIterator::new(BufReader::new(adapter_file));

    let read_padded = PaddedBytes::from_str::<NucMatrix>(seq, max_block_size);

    let mut local = Block::<true, false, false, true, true>::new(min_block_size, len, max_block_size);

    for line in iterator {
        let line = line?;
        let adapter = line.get_sequence();

        let mut cigar = Cigar::new(min_block_size, len);
        
        
        let adapter_padded = PaddedBytes::from_str::<NucMatrix>(adapter, max_block_size);
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

        let res_offset = AlignResult { score: res.score, query_idx: res.query_idx, reference_idx: res.reference_idx+offset };
        
        let result = Adapter::new(
            sequence.get_name().split_once(' ').unwrap().0.to_owned(),
            line.get_name().to_owned(),
            cigar.to_string(),
            res_offset
        );
        output.push(result);
    }
    Ok(())
}

fn remove_t(sequence: &mut Fasta) -> Option<AlignResult>{
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };
    let seq = sequence.get_sequence().to_owned();
    let seq_first = &seq[0..=300];

    let poly_t = 
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";


    let min_block_size = 256;
    let max_block_size = 512;
    
    let read_padded_first = PaddedBytes::from_str::<NucMatrix>(seq_first, max_block_size);
    let polyt_padded = PaddedBytes::from_str::<NucMatrix>(poly_t, max_block_size);
    
    
    let mut local = Block::<false, false, true, false, true>::new(poly_t.len(), 500, max_block_size);
    
    local.align(
        &polyt_padded, 
        &read_padded_first, 
        &NW1, 
        gaps, 
        min_block_size..=max_block_size, 
        0
    );



    let res_first = local.res();
    println!("{:?}", res_first);
    if res_first.score >= 11 {
        sequence.set_sequence(seq[res_first.reference_idx..].to_owned());
        return Some(res_first)
    } else {

        let seq_new = sequence.get_sequence().to_owned();
        let seq_reversed = reverse_complement(&seq_new);
        let seq_rev_first = &seq_reversed[..300];

        let read_padded_last = PaddedBytes::from_str::<NucMatrix>(seq_rev_first, max_block_size);


        local.align(
            &polyt_padded, 
            &read_padded_last, 
            &NW1, 
            gaps, 
            min_block_size..=max_block_size, 
            0
        );

        let res_last = local.res();

        if res_last.score >= 11 {
            sequence.set_sequence(seq_rev_first[res_last.reference_idx..].to_owned());
            return Some(res_last)
        }
    }
    None
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T', b'T' => 'A', b'C' => 'G', b'G' => 'C',
            b'a' => 't', b't' => 'a', b'c' => 'g', b'g' => 'c',
            _ => b as char,
        })
        .collect()
}


