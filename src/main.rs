// imports

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*imports][imports:1]]
#[macro_use]
extern crate clap;

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader};
use std::path::Path;

use clap::{App, AppSettings, Arg};
use trajectory_analysis::{analyze_frames, extract_frame};
// imports:1 ends here

// config

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*config][config:1]]
struct Config {
    inpfile  : Option<String>,
    extract  : bool,
    timestep : usize,
    outfile  : Option<String>,
    maxcols  : usize,
}

impl Config {
    fn new() -> Self {
        Config {
            inpfile  : None,
            extract  : false,
            timestep : 0,
            outfile  : None,
            maxcols  : 100,
        }
    }
}
// config:1 ends here

// fn parse_args

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*fn%20parse_args][fn parse_args:1]]
/// get file name from command line argument
fn parse_args() -> Result<Config, Box<dyn Error>> {
    let matches = App::new("trajectory-analysis")
        .version(crate_version!())
        .author(crate_authors!())
        .about("tools for lammps/reaxff reaction trajectory analysis")
        .arg(
            Arg::with_name("debug")
                .help("debug switch")
                .long("debug")
                .multiple(true)
                .short("d")
        )
        .arg(
            Arg::with_name("extract")
                .help("extract structure on a specific timestep")
                .long("extract")
                .short("e")
                .value_name("TIMESTEP")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("maxcols")
                .help("set max columns in fragment analysis report")
                .long("maxcols")
                .short("m")
                .value_name("NCOLUMNS")
                .takes_value(true)
        )
        .arg(
            Arg::with_name("input")
                .help("set input file name")
                .value_name("INPUT-FILE")
                .index(1)
                .required(true)
        )
        .arg(
            Arg::with_name("output")
                .help("set output file name")
                .value_name("OUTPUT-FILE")
                .index(2)
                .required(true)
        )
        .setting(AppSettings::ArgRequiredElseHelp)
        .get_matches();

    let mut config = Config::new();
    if let Some(r) = matches.value_of("extract") {
        config.extract = true;
        match r.parse::<usize>() {
            Ok(v) => {
                config.timestep = v;
            }
            Err(why) => {
                let msg = format!("timestep is expected to be integer, but found: {}", r);
                Err(msg)?;
            },
        }
    }

    if let Some(r) = matches.value_of("input") {
        config.inpfile = Some(r.to_string());
    } else {
        Err("input file is required!")?;
    }

    if let Some(r) = matches.value_of("output") {
        config.outfile = Some(r.to_string());
    } else {
        Err("output file is required!")?;
    }

    if let Some(r) = matches.value_of("maxcols") {
        config.maxcols = r.parse().unwrap();
    }

    Ok(config)
}
// fn parse_args:1 ends here

// fn main

// [[file:~/Workspace/Programming/structure-predication/trajectory-analysis/trajectory.note::*fn%20main][fn main:1]]
fn main() {
    let config = parse_args().unwrap();
    let inpfile = config.inpfile.unwrap();
    let outfile = config.outfile.unwrap();

    if config.extract {
        println!("extract structure on timestep {} from file: {}", config.timestep, inpfile);

        match extract_frame(&inpfile, config.timestep, &outfile) {
            Ok(r) => {
                println!("Done. CIF file saved as {}.", outfile);
            },
            Err(why) => {
                println!("Failed: {:}", why);
            },
        }
    } else {
        println!("fragment analysis of trajectory file: {:}", inpfile);
        match analyze_frames(&inpfile, &outfile, config.maxcols) {
            Ok(r) => {
                println!("Done. report file saved as {}", outfile);
            },
            Err(why) => {
                println!("Failed: {:}", why);
            },
        }
    }
}
// fn main:1 ends here
