open Core
open Bistro
open Bistro.Shell_dsl
open Bistro_utils

let img = [ docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" () ]

let input_dep p = dep (Workflow.input p)

let init ~dbname ~rootdir ~email ~famprefix ~ncbi_assembly_dir ~ncbi_taxonomy_dir ~custom_assembly_dir ~ncbi_assembly_4annot_dir =
  let r x = Filename.concat rootdir x in
  Workflow.shell ~descr:"pantagruel.init" [
    cmd "/pipeline/pantagruel_pipeline_init.sh" ~img [
      string dbname ;
      string rootdir ;
      string "/pantagruel" ;
      string email ;
      string famprefix ;
      string (Option.value ~default:(r"NCBI_assembly") ncbi_assembly_dir) ;
      string (Option.value ~default:(r"NCBI_taxonomy") ncbi_taxonomy_dir) ;
      string (Option.value ~default:(r"user_genome") custom_assembly_dir) ;
      string (Option.value ~default:(r"NCBI_assembly_annot") ncbi_assembly_4annot_dir) ;
    ]
  ]

let main ~outdir ~ncbi_assembly_dir ~ncbi_taxonomy_dir ~custom_assembly_dir ~ncbi_assembly_4annot_dir () =
  let w = init ~dbname:"ngono" ~rootdir:"." ~email:"philippe.veber@univ-lyon1.fr" ~famprefix:"NGONO" ~ncbi_assembly_dir ~ncbi_taxonomy_dir ~custom_assembly_dir ~ncbi_assembly_4annot_dir in
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
