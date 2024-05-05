# Fetch sequences from `tblout` nhmmer output

An automated way of getting the sequences from `tblout` nhmmer file. Under the hood, uses `esl-sfetch` and iterates over the `tblout` file.

## Usage

May change at any stage.

```console
Running extract_nhmmer_tblout
Usage: extract_nhmmer_tblout [OPTIONS] --esl-sfetch <esl-sfetch> <TBL> <FASTA>

Arguments:
  <TBL>    Path to the nhmmer tblout file.
  <FASTA>  Path to the fasta file used for nhmmer output.

Options:
  -e, --esl-sfetch <esl-sfetch>
          Path to esl-sfetch. If not installed, it's part of HMMER.
  -v, --e-value-threshold <e_value_threshold>
          E-value threshold for hits to extract. [default: 0.00001]
  -h, --help
          Print help
  -V, --version
          Print version
```

## Requirements

You'll need the `easel` part of HMMER, and point to the executable (`-v /path/to/esl-sfetch`).
