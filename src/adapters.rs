use std::{cmp::min, fs::File, io::BufReader};

use anyhow::{anyhow, bail, Ok, Result};
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

//Aligns adapter sequence against ref sequence
pub fn align_adapters(fasta: &mut Fasta, cmds: &AdapterArgs, output: &mut Vec<Adapter>) -> Result<()> {
    
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };

    let mut offset = 0;

    if cmds.remove_t{
        let t_res = remove_t(fasta, cmds.score_cutoff)?;
        offset = match t_res {
            Some(t_res) => t_res.reference_idx,
            None => 0
        };
    }
    
    let seq = fasta.get_sequence();
    let len = seq.len();

    let min_block_size: usize = percent_len(cmds.largest_adapter_len, 1.0);
    let max_block_size: usize = 4*min_block_size;

    let adapter_file = File::open(&cmds.adapters_file_path)?;

    let iterator = FastaIterator::new(BufReader::new(adapter_file));

    let read_padded = PaddedBytes::from_str::<NucMatrix>(seq, max_block_size);

    let mut local = Block::<true, false, false, true, true>::new(min_block_size, len, max_block_size);

    for line in iterator {
        let line = line?;
        
        let adapter = line.get_sequence();

        if cmds.largest_adapter_len < adapter.len(){
            return Err(anyhow!("Largest adapter length not large enough"));
        }

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
        if res.score >= cmds.score_cutoff {
            local.trace().cigar_eq(&adapter_padded, &read_padded, res.query_idx, res.reference_idx, &mut cigar);

            let res_offset = AlignResult { score: res.score, query_idx: res.query_idx, reference_idx: res.reference_idx+offset };
            
            let result = Adapter::new(
                fasta.get_name().split_once(' ').unwrap().0.to_owned(),
                line.get_name().to_owned(),
                cigar.to_string(),
                res_offset
            );
            output.push(result);
        }
    }
    Ok(())
}


//Removes poly t at start or poly a at end and everything before that
fn remove_t(fasta: &mut Fasta, smallest: i32) -> Result<Option<AlignResult>>{
    let gaps = Gaps {
        open: -2,
        extend: -1,
    };
    let seq = fasta.get_sequence().to_owned();
    let seq_first = &seq[0..min(300, seq.len())];

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

    if res_first.score >= smallest {
        fasta.set_sequence(seq[res_first.reference_idx..].to_owned());
        return Ok(Some(res_first))
    
    } else {

        let seq_new = fasta.get_sequence().to_owned();
        let seq_reversed = reverse_complement(&seq_new)?;
        let seq_rev_first = &seq_reversed[..min(seq_new.len(),300)];

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

        if res_last.score >= smallest {
            fasta.set_sequence(seq_rev_first[res_last.reference_idx..].to_owned());
            return Ok(Some(res_last))
        }
    }
    Ok(None)
}

fn reverse_complement(seq: &str) -> Result<String> {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => Ok('T'), 
            b'T' => Ok('A'), 
            b'C' => Ok('G'), 
            b'G' => Ok('C'),
            b'a' => Ok('t'), 
            b't' => Ok('a'), 
            b'c' => Ok('g'), 
            b'g' => Ok('c'),
            _ => bail!("Not a base pair - check format"),
        })
        .collect()
}


