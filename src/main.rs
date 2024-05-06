use std::{
    ffi::OsStr,
    fs::File,
    path::{Path, PathBuf},
    process::Stdio,
};

use anyhow::{Context, Result};
use clap::{arg, command, value_parser, Arg};
use hmm_tblout::Reader;
use std::process::Command as Cmd;
use tempfile::tempdir;

fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename).extension().and_then(OsStr::to_str)
}

fn main() -> Result<()> {
    eprintln!("Running extract_nhmmer_tblout");
    // set up the app
    let matches = command!()
        .arg_required_else_help(true)
        .arg(
            arg!(<TBL> "Path to the nhmmer tblout file.")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            arg!(<FASTA> "Path to the fasta file used for nhmmer output.")
                .required(true)
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("esl-sfetch")
                .short('e')
                .long("esl-sfetch")
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .help("Path to esl-sfetch. If not installed, it's part of HMMER."),
        )
        .arg(
            Arg::new("e_value_threshold")
                .short('v')
                .long("e-value-threshold")
                .value_parser(value_parser!(f32))
                .required(false)
                .default_value("0.00001")
                .help("E-value threshold for hits to extract."),
        )
        .get_matches();

    // get the matches
    let tbl = matches
        .get_one::<PathBuf>("TBL")
        .expect("tbl is required")
        .clone();

    let fasta = matches
        .get_one::<PathBuf>("FASTA")
        .expect("fasta is required")
        .clone();

    let esl_sfetch = matches
        .get_one::<PathBuf>("esl-sfetch")
        .expect("esl-sfetch is required")
        .clone();

    let e_value_threshold = *matches
        .get_one::<f32>("e_value_threshold")
        .expect("defaulted by clap");

    // copy the fasta to a temporary directory
    let tmpdir = tempdir().context("Could not create tempdir")?;

    // check if the fasta is gzipped
    // if it is, use gunzip -c to copy to tmpdir
    // else just copy over
    let fasta_is_gzipped =
        get_extension_from_filename(fasta.to_str().context("Could not convert path to string")?)
            == Some("gz");
    let new_fasta_path = if fasta_is_gzipped {
        eprintln!("Input fasta is gzipped, unzipping...");

        let fasta_file_name = fasta
            .clone()
            .file_stem()
            .context("Could not get file stem")?
            .to_os_string();

        let fasta_file = File::create(tmpdir.path().join(&fasta_file_name))
            .context("Could not create fasta file")?;
        let stdio = Stdio::from(fasta_file);
        let copy_via_gzip = Cmd::new("gunzip")
            .arg("-c")
            .arg(fasta.clone())
            .stdout(stdio)
            .spawn()?;
        copy_via_gzip.wait_with_output()?;

        fasta_file_name
    } else {
        eprintln!("Input fasta is not gzipped, copying...");
        let copy_over = Cmd::new("cp")
            .arg(fasta.clone())
            .arg(tmpdir.path())
            .spawn()?;
        copy_over.wait_with_output()?;

        let fasta_file_name = fasta.clone();
        fasta_file_name.into_os_string()
    };

    // index the fasta
    let new_fasta_location = tmpdir.path().join(new_fasta_path.clone());
    eprintln!("New fasta location: {:?}", new_fasta_location);
    eprintln!("Indexing fasta");
    let index_fasta = Cmd::new(esl_sfetch.clone())
        .arg("--index")
        .arg(new_fasta_location.clone())
        .spawn()?;
    index_fasta.wait_with_output()?;

    // that takes a moment.
    // now we can open the tbl, iterate over, and extract the sequences
    let mut reader = Reader::from_path(tbl)?;

    eprintln!("Iterating over tblout");
    for record in reader.records() {
        let r = record?;

        // not interested in low value hits
        if r.e_value() > e_value_threshold {
            continue;
        }

        let target_name = r.target_name();
        let ali_from_to = format!("{}..{}", r.ali_from(), r.ali_to());

        let extract_sequences = Cmd::new(esl_sfetch.clone())
            .arg("-c")
            .arg(ali_from_to)
            .arg(new_fasta_location.clone())
            .arg(target_name)
            .output()?;

        // TODO: maybe modify the fasta header to include extra information
        let out = String::from_utf8(extract_sequences.stdout)?;

        print!("{}", out);
    }

    Ok(())
}
