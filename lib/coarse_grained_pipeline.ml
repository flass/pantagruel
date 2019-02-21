open Core
open Bistro
open Bistro.Shell_dsl
open Bistro_utils

let img = [ docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" () ]

let input_dep p = dep (Workflow.input p)

let env xs =
  list
    ~sep:"\n"
    (fun (k, v) -> seq ~sep:"" [ string (sprintf "export %s=" k) ; v ])
    xs

let environment
    ?(ptgrepo = string "/pantagruel")
    ?(email = string "undisclosed")
    ?(famprefix = string "PANTAG")
    ?(customassemb = dest // "user_genome")
    ?(ncbiass = dest // "NCBI/Assembly")
    ?(refass = dest // "NCBI/Assembly")
    ?(ncbitax = dest // "NCBI/Taxonomy")
    ?(chaintype = string "fullgenetree")
    ?(pseudocoremingenomes = string "")
    ?(coreseqtype = string "cds")
    ?(poplgthresh = string "default")
    ?(popbsthresh = string "default")
    ?(hpcremoteptgroot = string "none")
    ~ptgroot ~ptgdbname
    ()
  =
  env [
    "ptgrepo", ptgrepo ;
    "email", email ;
    "famprefix", famprefix ;
    "customassemb", customassemb ;
    "ncbiass", ncbiass ;
    "refass", refass ;
    "ncbitax", ncbitax ;
    "chaintype", chaintype ;
    "pseudocoremingenomes", pseudocoremingenomes ;
    "coreseqtype", coreseqtype ;
    "poplgthresh", poplgthresh ;
    "popbsthresh", popbsthresh ;
    "hpcremoteptgroot", hpcremoteptgroot ;
    "ptgroot", ptgroot ;
    "ptgdbname", ptgdbname ;
  ]

let script xs = seq ~sep:"\n" (List.map ~f:(seq ~sep:" ") xs)

let init ~ptgdbname ~ptgroot ?email ?famprefix ~ncbi_assembly_dir:_ ~ncbi_taxonomy_dir:_ ~custom_assembly_dir:_ ~ncbi_assembly_4annot_dir:_ () =
  let env =
    environment
      ~ptgdbname ~ptgroot ?email ?famprefix ()
  in
  let script = script [
      [ string "source" ; file_dump env ] ;
      [ string "/pantagruel/scripts/pipeline/pantagruel_pipeline_init.sh" ] ;
    ]
  in
  Workflow.shell ~descr:"pantagruel.init" [
    cmd "bash" ~img [ file_dump script ]
  ]

let main ~outdir ~ncbi_assembly_dir ~ncbi_taxonomy_dir ~custom_assembly_dir ~ncbi_assembly_4annot_dir () =
  let w =
    init
      ~ptgdbname:(string "ngono")
      ~ptgroot:dest
      ~ncbi_assembly_dir
      ~ncbi_taxonomy_dir
      ~custom_assembly_dir
      ~ncbi_assembly_4annot_dir
      ()
  in
  Repo.build_main ~outdir ~np:4 ~mem:(`GB 4) Repo.[ item ["delme"] w ]

let command =
  let open Command.Let_syntax in
  Command.basic
    ~summary:"coarse-grained pipeline"
    [%map_open
      let outdir = flag "--outdir" (required string) ~doc:"PATH Destination directory."
      and ncbi_assembly_4annot_dir = flag "--ncbi-assembly-4annot-dir" (optional string) ~doc:"PATH ncbi_assembly_4annot_dir"
      and ncbi_taxonomy_dir = flag "--ncbi-taxonomy-dir" (optional string) ~doc:"PATH ncbi_taxonomy_dir"
      and custom_assembly_dir = flag "--custom-assembly-dir" (optional string) ~doc:"PATH custom_assembly_dir"
      and ncbi_assembly_dir = flag "--ncbi-assembly-dir" (optional string) ~doc:"PATH ncbi_assembly_dir"
      in
      main ~outdir ~ncbi_assembly_4annot_dir ~ncbi_taxonomy_dir ~custom_assembly_dir ~ncbi_assembly_dir
    ]
