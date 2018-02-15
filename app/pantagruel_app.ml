open Core
open Pantagruel
open Bistro_utils

let run ~outdir ~taxid () =
  let pipeline = Pipeline.make ~taxid in
  let repo = Pipeline.repo pipeline in
  (* Lwt_main.run (Entrez.assembly_request ~taxid) ; *)
  Repo.build ~outdir repo

let command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"Pantagruel"
    [%map_open
      let outdir =
        flag "--outdir" (required string) ~doc:"PATH Destination directory."
      and taxid =
        flag "--taxid"  (required int) ~doc:"INTEGER NCBI taxid" in
      fun () ->
        run ~outdir ~taxid ()
    ]

let () = Command.run ~version:"dev" command
