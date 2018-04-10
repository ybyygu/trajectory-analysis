// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::68b8f3aa-b3f8-43c0-8b4d-c3165b146535][68b8f3aa-b3f8-43c0-8b4d-c3165b146535]]
extern crate reaxa;
#[macro_use]
extern crate clap;

use std::fs::File;
use std::error::Error;
use std::io::{self, BufReader};
use std::io::prelude::*;
use std::collections::HashMap;
use std::path::Path;

use clap::{App, Arg, AppSettings};
use reaxa::{extract_frame, analyze_frames};
// 68b8f3aa-b3f8-43c0-8b4d-c3165b146535 ends here

// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::25390ea6-8bf3-4083-a422-f9e628299efc][25390ea6-8bf3-4083-a422-f9e628299efc]]
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
// 25390ea6-8bf3-4083-a422-f9e628299efc ends here

// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::b8ea57f0-b549-4fa0-ac1a-abf83009009e][b8ea57f0-b549-4fa0-ac1a-abf83009009e]]
/// get file name from command line argument
fn parse_args() -> Result<Config, Box<Error>> {
    let matches = App::new("reaxa")
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
// b8ea57f0-b549-4fa0-ac1a-abf83009009e ends here

// [[file:~/Workspace/Programming/structure-predication/reaction-analysis/reaction-analysis.note::0176a8ca-12d4-4334-9fa5-6093ceddb854][0176a8ca-12d4-4334-9fa5-6093ceddb854]]
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
// 0176a8ca-12d4-4334-9fa5-6093ceddb854 ends here
