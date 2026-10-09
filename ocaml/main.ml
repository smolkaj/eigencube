open Base
open Stdio
open Eigencube

let print_usage () =
  printf "Usage: dune exec eigencube -- [seed]\n";
  printf
    "  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)\n"

let () =
  let args = Sys.get_argv () in
  if Array.length args > 1 then
    match args.(1) with
    | "--help" | "-h" -> print_usage ()
    | s ->
    match Int.of_string_opt s with
    | Some seed ->
      log "Solving scrambled cube (seed=%d)..." seed;
      ignore (solve (shuffle (solved_cube ()) 100_000 seed))
    | None ->
      printf "Error: unrecognized option or invalid integer seed '%s'.\n\n%!" s;
      print_usage ();
      Stdlib.exit 1
  else begin
    log "Solving scrambled cube (seed=42)...";
    ignore (solve (shuffle (solved_cube ()) 100_000 42))
  end
