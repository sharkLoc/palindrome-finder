use std::{
    cmp::{max, min},
    mem,
};

use crate::{
    command_line::WfaArgs,
    fasta_parsing::Fasta,
    output::PalindromeData,
};

use anyhow::{bail, ensure, Ok, Result};

const SIZE: usize = 1000;

//Use WFA algorithm to find palindromes
pub fn wfa_palins(
    fasta: Fasta,
    output: &mut Vec<PalindromeData>,
    wfa_args: &WfaArgs,
) -> Result<()> {

    //Converts the sequence to bytes, where A = !T, C = !G
    let mut seq_clone = fasta.sequence.clone();
    let bytes_seq = unsafe { seq_clone.as_bytes_mut() };
    sequence_to_bytes(bytes_seq)?;

    let seq = fasta.sequence;

    let len = seq.len();
    let mut index = 0;

    let mut wf = vec![0; max(SIZE, wfa_args.gap_len + 2)];
    let mut wf_next = vec![0; max(SIZE, wfa_args.gap_len + 2)];

    let first_wave = vec![0; wfa_args.gap_len + 2];

    while index <= len {
        let mut edit_dist = 0;
        let mut wf_len = wfa_args.gap_len + 1;

        let mut max_index = 0;
        let mut max_score = 0.0;

        //Reset first wave to 0s
        wf[..=wf_len].copy_from_slice(&first_wave);

        'outer: while (edit_dist as f32) / (wf[max_index] as f32 + 0.001)
            <= wfa_args.mismatch_proportion
        {
            let mut max_wf_score = 0.0;
            for i in 0..wf_len {

                //Extend wave along the matches
                let (mut x, mut y) = get_xy(wf_len, i, wf[i], wfa_args.gap_len);
                x += index;
                let counter = extend_wave(x, y, index, bytes_seq)?; //Using bytes_seq for fast complement checks

                wf[i] += counter as usize;
                x += counter as usize;
                y += counter as usize;

                let score = calculate_score(x, y, edit_dist, wfa_args);
                max_wf_score = f32::max(max_wf_score, score);

                if wf[i] > wf[max_index] {
                    max_index = i;
                }

                if x == len || y == index {
                    break 'outer;
                }
            }

            max_score = f32::max(max_score, max_wf_score);

            //X-drop pruning
            if max_wf_score < max_score - f32::max(wfa_args.x_drop, (wf[max_index] as f32) * 0.1) {
                break;
            }

            next_wave(&mut wf, &mut wf_next, wf_len);
            max_index += 1;
            edit_dist += 1;
            wf_len += 2;
        }

        if wf[max_index] == 0 {
            index += 1;
            continue;
        }

        let mut increment = 1;

        let (x, y) = get_xy(wf_len, max_index, wf[max_index], wfa_args.gap_len);
        let palin = &seq[index - y..index + x];
        let gap = y - wf[max_index];

        if x >= wfa_args.min_length {
            let palin = PalindromeData::new(
                (index - y) as u32,
                (index + x - 1) as u32,
                x as u32,
                gap as u32,
                (x + y) as u32,
                edit_dist,
                fasta.name.to_owned(),
                palin.to_owned(),
            );
            output.push(palin);
            increment = x;
        }
        index += increment
    }
    Ok(())
}

fn calculate_score(x: usize, y: usize, d: u32, args: &WfaArgs) -> f32 {
    (x + y) as f32 * (args.match_bonus / 2.0)
        - (d as f32) * (args.match_bonus - (-args.mismatch_penalty))
}

fn sequence_to_bytes(seq: &mut [u8]) -> Result<()> {
    for i in seq.iter_mut() {
        *i = match i {
            65 | 97 => 2,    // A, a
            84 | 116 => 253, // T, t
            67 | 99 => 3,    // C, c
            71 | 103 => 252, // G, g
            _ => bail!("Not a base pair - check format"),
        };
    }
    Ok(())
}

//Extends wave across matches
fn extend_wave(mut x: usize, mut y: usize, index: usize, seq: &[u8]) -> Result<u32> {
    let len = seq.len();
    let mut count = 8;
    let mut counter = 0;

    while x < len && y < index && count >= 8 {
        let len1 = min(len - x, 8);
        let len2 = min(index - y, 8);
        count = count_matching(&seq[x..x + len1], &seq[index - y - len2..index - y])?;
        x += count as usize;
        y += count as usize;
        counter += count;
    }
    Ok(counter)
}

//Counts the matching sequences with bit manipulation
fn count_matching(seq1: &[u8], seq2: &[u8]) -> Result<u32> {
    ensure!(
        seq1.len() <= 8,
        "Sequence length too long when processing bits"
    );
    ensure!(
        seq2.len() <= 8,
        "Sequence length too long when processing bits"
    );
    let mut buf1 = [0; 8];
    let mut buf2 = [0; 8];
    buf1[..seq1.len()].copy_from_slice(seq1);
    buf2[8 - seq2.len()..].copy_from_slice(seq2);
    let num1 = u64::from_le_bytes(buf1);
    let num2 = !u64::from_be_bytes(buf2);
    let diff = num1 ^ num2;
    Ok(diff.trailing_zeros() / 8)
}

//Branches out to next wave
fn next_wave(wf: &mut Vec<usize>, wf_next: &mut Vec<usize>, wf_len: usize) {
    for i in 0..wf_len {
        if i == 0 {
            wf_next[i] = wf[i];
            wf_next[i + 1] = max(wf[i] + 1, wf[i + 1]);
        } else if i != wf_len - 1 {
            if i + 1 >= wf_next.len() {
                wf_next.extend([0; SIZE]);
                wf.extend([0; SIZE])
            }
            wf_next[i + 1] = max(wf[i] + 1, max(wf[i - 1], wf[i + 1]));
        } else {
            if i + 1 >= wf_next.len() || i + 2 >= wf_next.len() {
                wf_next.extend([0; SIZE]);
                wf.extend([0; SIZE])
            }
            wf_next[i + 2] = wf[i];
            wf_next[i + 1] = max(wf[i - 1], wf[i] + 1);
        }

        if wf_len == 1 {
            wf_next[i + 2] = wf[i];
        }
    }
    mem::swap(wf, wf_next);
}

//Computes XY coordinates of a wave
fn get_xy(wf_len: usize, index: usize, length: usize, gap_size: usize) -> (usize, usize) {
    let offset = wf_len as i32 - ((wf_len - (gap_size + 1)) / 2) as i32 - (index as i32) - 1;
    if offset >= 0 {
        (length, length + offset as usize)
    } else {
        (length + offset.unsigned_abs() as usize, length)
    }
}