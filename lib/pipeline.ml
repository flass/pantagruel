open Core
open Bistro.Std
open Bistro_bioinfo.Std
open Bistro.EDSL
open Bistro_utils

let tools_env = docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" ()

let collect_all_proteins assembly_folder =
  workflow ~descr:"collect_all_proteins" [
    cmd "bash" ~env:tools_env [
      file_dump (string Scripts.collect_all_proteins) ;
      dep assembly_folder ;
      dest ;
    ] ;
  ]

let dereplicate_fasta (fa : fasta workflow) : fasta workflow =
  workflow ~descr:"dereplicate_fasta" [
    cmd "python" ~env:tools_env [
      file_dump (string Scripts.dereplicate_fasta) ;
      dep fa ;
      dest ;
    ] ;
  ]

let split_mmseqs_clustdb_fasta fa =
  workflow ~descr:"split_mmseqs_clustdb_fasta" [
    cmd "python" ~env:tools_env [
      file_dump (string Scripts.split_mmseqs_clustdb_fasta) ;
      dep fa ;
      string "ABCDEP" ;
      dest ;
    ]
  ]

let make assembly_folder_path =
  let assembly_folder = input assembly_folder_path in
  let all_unique_proteins =
    collect_all_proteins assembly_folder
    |> dereplicate_fasta
  in
  let protein_families =
    let open Mmseqs2 in
    let db = createdb all_unique_proteins in
    createseqfiledb db (cluster db) in
  object
    method output = split_mmseqs_clustdb_fasta protein_families
  end


let repo p =
  Repo.[
    item ["output"] p#output ;
  ]
