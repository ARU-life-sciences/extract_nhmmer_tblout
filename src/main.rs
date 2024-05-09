use std::{
    ffi::OsStr,
    fs::File,
    io,
    path::{Path, PathBuf},
    process::Stdio,
};

use anyhow::{Context, Result};
use clap::{arg, command, value_parser, Arg};
use fasta::record::Definition;
use hmm_tblout::Reader;
use noodles_fasta as fasta;
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
            arg!([FASTA] "Path to the fasta file used for nhmmer output. If not specified, the target file from the tblout file is used (this probably only works when that file path is absolute).")
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
        .arg(
            Arg::new("species_id")
                .short('s')
                .long("species-id")
                .value_parser(value_parser!(String))
                .required(false)
                .default_value("")
                .help("Species ID to add to the start of the header. Useful for downstream processing."),
        )
        .get_matches();

    // get the matches
    let tbl = matches
        .get_one::<PathBuf>("TBL")
        .expect("tbl is required")
        .clone();

    let fasta_match = matches.get_one::<PathBuf>("FASTA").cloned();

    let esl_sfetch = matches
        .get_one::<PathBuf>("esl-sfetch")
        .expect("esl-sfetch is required")
        .clone();

    let e_value_threshold = *matches
        .get_one::<f32>("e_value_threshold")
        .expect("defaulted by clap");

    let mut species_id = matches
        .get_one::<String>("species_id")
        .expect("defaulted by clap")
        .clone();

    // copy the fasta to a temporary directory
    let tmpdir = tempdir().context("Could not create tempdir")?;

    // read the tblout to ge the metadata
    let mut reader = Reader::from_path(tbl)?;
    let meta = reader.meta();
    let target_file = meta.target_file();

    let fasta = match fasta_match {
        Some(f) => f,
        None => target_file,
    };

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

    eprintln!("Iterating over tblout");
    for record in reader.records() {
        let r = record?;

        // not interested in low value hits
        if r.e_value().unwrap() > e_value_threshold {
            continue;
        }

        let target_name = r.target_name();
        let ali_from_to = format!("{}..{}", r.ali_from().unwrap(), r.ali_to().unwrap());

        let extract_sequences = Cmd::new(esl_sfetch.clone())
            .arg("-c")
            .arg(ali_from_to)
            .arg(new_fasta_location.clone())
            .arg(target_name)
            .output()?;

        // parse the fasta properly and edit the header.
        let mut parsed_fasta = fasta::reader::Reader::new(&extract_sequences.stdout[..]);
        let stdout = io::stdout().lock();
        let mut writer = fasta::Writer::new(stdout);

        for record in parsed_fasta.records() {
            let r = record?;

            let append_name = std::str::from_utf8(r.name())?;
            species_id += append_name;

            let def = Definition::new(species_id.as_bytes(), r.description().map(|e| e.to_vec()));

            let new_record = fasta::Record::new(def, r.sequence().to_owned());
            writer.write_record(&new_record)?;
        }
    }

    Ok(())
}
