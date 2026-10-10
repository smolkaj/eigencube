open Base
open Stdio
open Eigencube

let print_usage () =
  printf "Usage: dune exec eigencube -- [seed]\n";
  printf
    "  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)\n"

let () =
  match Sys.get_argv () with
  | [| _ |] ->
    log "Solving scrambled cube (seed=42)...";
    ignore (solve (shuffle ~iters:100_000 ~seed:42 solved_cube))
  | [| _; "--help" | "-h" |] -> print_usage ()
  | [| _; s |] -> (
    match Int.of_string_opt s with
    | Some seed ->
      log "Solving scrambled cube (seed=%d)..." seed;
      ignore (solve (shuffle ~iters:100_000 ~seed solved_cube))
    | None ->
      printf "Error: unrecognized option or invalid integer seed '%s'.\n\n%!" s;
      print_usage ();
      Stdlib.exit 1
  )
  | _ ->
    printf "Error: unexpected extra arguments.\n\n%!";
    print_usage ();
    Stdlib.exit 1
