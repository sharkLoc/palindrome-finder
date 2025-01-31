use std::fmt::Display;

use clap::{command, ArgGroup, Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about)]
pub struct PalinArgs {
    ///Decide which algorithm should be used
    #[clap(subcommand)]
    pub mode: AlgorithmType,
}

#[derive(Debug, Subcommand)]
pub enum AlgorithmType {
    ///Use fast WFA algorithm, allows mismatches and indels and uses more complex pruning
    Wfa(WfaArgs),
    ///Use fixed-mismatches algorithm, only allows fixed number mismatches and no indels
    ExactMatch(FixedArgs), 
    ///Script for aligning adapter sequences, uses block-align library
    Adapters(AdapterArgs)
}

#[derive(Debug, Args)]
#[command(group = ArgGroup::new("file_type")
    .required(true)
    .args(&["fa", "fgz", "fq", "fqgz"]))]
pub struct AdapterArgs{
    #[arg(short, long = "input", required = true)]
    ///Input file path
    pub input_file: String,

    /// Indicates the input file should be read in FASTA format
    #[arg(long)]
    pub fa: bool,

    /// Indicates the input file should be read in compressed FASTA gzip format
    #[arg(long)]
    pub fgz: bool,

    ///Indicates the input file should be read in FASTQ format
    #[arg(long)]
    pub fq: bool,

    ///Indicates the input file should be read in compressed FASTQ gzip format
    #[arg(long)]
    pub fqgz: bool,

    #[arg(short, long = "output")]
    ///Output file path.
    pub output_file: String,

    ///The file path for the list of adapter seqeuences. Must be fasta format
    #[arg(short, long)]
    pub adapters_file_path: String,
    
    ///The length of the longest adapter in the file
    #[arg(short, long)]
    pub longest_adapter: usize,

    ///The smallest alignment score for an adapter sequence to be outputted
    #[arg(short = 'c', long)]
    pub score_cutoff: i32,

    ///Enables removing  poly-t sequences at the start and all the sequences before.
    ///This also removes any poly-a sequences at the end if no poly-t are found
    #[arg(long)]
    pub remove_t: bool

}

#[derive(Debug, Args)]
#[command(group = ArgGroup::new("file_type")
    .required(true)
    .args(&["fa", "fgz", "fq", "fqgz"]))]
pub struct FixedArgs{
    #[arg(short, long = "input", required = true)]
    ///Input file path
    pub input_file: String,

    /// Indicates the input file should be read in FASTA format
    #[arg(long)]
    pub fa: bool,

    /// Indicates the input file should be read in compressed FASTA gzip format
    #[arg(long)]
    pub fgz: bool,

    ///Indicates the input file should be read in FASTQ format
    #[arg(long)]
    pub fq: bool,

    ///Indicates the input file should be read in compressed FASTQ gzip format
    #[arg(long)]
    pub fqgz: bool,

    #[arg(short, long = "output")]
    ///Output file path.
    pub output_file: String,

    #[arg(short = 'l', long, default_value_t = 10)]
    ///Minimum palindrome arm length
    pub len: usize,

    #[arg(short, long = "gap", default_value_t = 3)]
    ///Maximum gap length in a palindrome
    pub gap_len: usize,
    ///Max number of mismatches
    #[arg(short = 'm', long = "mismatches", default_value_t = 4)]
    pub mismatches: u32,
}

#[derive(Debug, Args)]
#[command(group = ArgGroup::new("file_type")
    .required(true)
    .args(&["fa", "fgz", "fq", "fqgz"]))]
pub struct WfaArgs {
    #[arg(short, long = "input", required = true)]
    ///Input file path
    pub input_file: String,

    /// Indicates the input file should be read in FASTA format
    #[arg(long)]
    pub fa: bool,

    /// Indicates the input file should be read in compressed FASTA gzip format
    #[arg(long)]
    pub fgz: bool,

    ///Indicates the input file should be read in FASTQ format
    #[arg(long)]
    pub fq: bool,

    ///Indicates the input file should be read in compressed FASTQ gzip format
    #[arg(long)]
    pub fqgz: bool,

    #[arg(short, long = "output")]
    ///Output file path.
    pub output_file: String,

    #[arg(short = 'l', long, default_value_t = 10)]
    ///Minimum palindrome arm length
    pub min_length: usize,

    #[arg(short, long = "gap", default_value_t = 3)]
    ///Maximum gap length in a palindrome
    pub gap_len: usize,
    
    ///Bonus for matches in scoring, must be positive
    #[arg(short = 'b', long, default_value_t = 1.0)]
    pub match_bonus: f32,

    ///Penalty for mismatches in scoring, must be positive (since value is subtracted)
    #[arg(short = 'p', long, default_value_t = 4.0)]
    pub mismatch_penalty: f32,

    ///Maximum score drop allowed before pruning
    #[arg(short, long, default_value_t = 20.0)]
    pub x_drop: f32,

    ///Max percentage of mismatches allowed in a palindrome, must be between 0 and 1
    #[arg(short = 'm', long, default_value_t = 0.05)]
    pub mismatch_proportion: f32,
}

impl Display for PalinArgs{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.mode
        )
    }
}

impl Display for AlgorithmType{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            AlgorithmType::Wfa(cmds) => 
                write!(
                    f,
                    "Min length: {}\nMax gap length: {}\nMatch bonus: {}\nMismatch penalty: {}\nX-drop: {}\nMax mismatch proportion: {}",
                    cmds.min_length, cmds.gap_len, cmds.match_bonus, cmds.mismatch_penalty, cmds.x_drop, cmds.mismatch_proportion
            ),
            AlgorithmType::ExactMatch(cmds) => 
                write!(
                    f,
                    "Min length: {}\nMax gap length: {}\nMismatches allowed: {}",
                    cmds.len, cmds.gap_len, cmds.mismatches
            ),
            AlgorithmType::Adapters(_cmds) => Ok(())
        }
        
    }
}

impl AlgorithmType {
    pub fn input_file(&self) -> &str {
        match self {
            AlgorithmType::Wfa(cmds) => &cmds.input_file,
            AlgorithmType::ExactMatch(cmds) => &cmds.input_file,
            AlgorithmType::Adapters(cmds) => &cmds.input_file,
        }
    }

    pub fn is_fa(&self) -> bool {
        match self {
            AlgorithmType::Wfa(cmds) => cmds.fa,
            AlgorithmType::ExactMatch(cmds) => cmds.fa,
            AlgorithmType::Adapters(cmds) => cmds.fa,
        }
    }

    pub fn is_fgz(&self) -> bool {
        match self {
            AlgorithmType::Wfa(cmds) => cmds.fgz,
            AlgorithmType::ExactMatch(cmds) => cmds.fgz,
            AlgorithmType::Adapters(cmds) => cmds.fgz,
        }
    }

    pub fn is_fq(&self) -> bool {
        match self {
            AlgorithmType::Wfa(cmds) => cmds.fq,
            AlgorithmType::ExactMatch(cmds) => cmds.fq,
            AlgorithmType::Adapters(cmds) => cmds.fq,
        }
    }

    pub fn is_fqgz(&self) -> bool {
        match self {
            AlgorithmType::Wfa(cmds) => cmds.fqgz,
            AlgorithmType::ExactMatch(cmds) => cmds.fqgz,
            AlgorithmType::Adapters(cmds) => cmds.fqgz,
        }
    }
    pub fn output_file(&self) -> &str {
        match self {
            AlgorithmType::Wfa(cmds) => &cmds.output_file,
            AlgorithmType::ExactMatch(cmds) => &cmds.output_file,
            AlgorithmType::Adapters(cmds) => &cmds.output_file,
        }
    }
}

