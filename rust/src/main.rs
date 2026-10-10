use eigencube::{log, shuffle, solve, SOLVED_CUBE};
use std::env;

fn print_usage() {
  println!("Usage: cargo run --release -- [seed]");
  println!("  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)");
}

fn main() {
  let args: Vec<String> = env::args().collect();
  match args.len() {
    1 => {
      log("Solving scrambled cube (seed=42)...");
      let scrambled = shuffle(100_000, 42, SOLVED_CUBE);
      solve(scrambled);
    }
    2 if args[1] == "--help" || args[1] == "-h" => {
      print_usage();
    }
    2 => match args[1].parse::<u64>() {
      Ok(seed) => {
        log(&format!("Solving scrambled cube (seed={})...", seed));
        let scrambled = shuffle(100_000, seed, SOLVED_CUBE);
        solve(scrambled);
      }
      Err(_) => {
        eprintln!("Error: unrecognized option or invalid integer seed '{}'.\n", args[1]);
        print_usage();
        std::process::exit(1);
      }
    },
    _ => {
      eprintln!("Error: unexpected extra arguments.\n");
      print_usage();
      std::process::exit(1);
    }
  }
}
