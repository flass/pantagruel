open Core
open Bistro.Std
open Bistro_bioinfo.Std
open Bistro.EDSL
open Bistro_utils

let tools_env = docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" ()

let assembly_list_of_folder (folder : _ directory workflow) =
  workflow ~descr:"assembly_list_of_folder" [
    cmd "ls" ~stdout:dest [
      string "-d" ;
      dep folder // "GC[AF]_*" ;
    ] ;
  ]

let collect_all_proteins (assembly_list : #text_file workflow)  =
  workflow ~descr:"collect_all_proteins" [
    cmd "bash" ~env:tools_env [
      file_dump (string Scripts.collect_all_proteins) ;
      dep assembly_list ;
      dest ;
    ] ;
  ]

let generate_assembly_stats assembly_folder =
  workflow ~descr:"generate_assembly_stats" [
    cmd "bash" ~env:tools_env [
      file_dump (string Scripts.generate_assembly_stats) ;
      dep assembly_folder ;
      dest ;
    ] ;
  ]

let extract_metadata_from_gbff ~assembly_list ~assembly_stats =
  workflow ~descr:"extract_metadata_from_gbff" [
    cmd "python" [
      file_dump (string Scripts.extract_metadata_from_gbff) ;
      opt "--assembly_folder_list" dep assembly_list ;
      opt "--add_assembly_info_dir" dep assembly_stats ;
      opt "--default_species_name" (string % quote ~using:'\'') "unclassified organism" ;
      opt "--output" ident dest ;
    ]
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

let stage1 assembly_folder_path =
  let assembly_folder : _ directory workflow = input assembly_folder_path in
  let assembly_list = assembly_list_of_folder assembly_folder in
  let assembly_stats = generate_assembly_stats assembly_folder in
  let assembly_metadata =
    extract_metadata_from_gbff ~assembly_list ~assembly_stats in
  let all_unique_proteins =
    collect_all_proteins assembly_list
    |> dereplicate_fasta
  in
  let protein_families =
    let open Mmseqs2 in
    let db = createdb all_unique_proteins in
    createseqfiledb db (cluster db) in
  object
    method assembly_list = assembly_list
    method assembly_stats = assembly_stats
    method assembly_metadata = assembly_metadata
    method protein_families = split_mmseqs_clustdb_fasta protein_families
  end


let stage2_pipeline fa =
  let alignment = Clustalo.clustalo fa in
  object
    method alignment = alignment
  end
  
let stage2 family_folder =
  let family_fastas =
    Sys.readdir family_folder
    |> Array.to_list
    |> List.map ~f:input
  in
  object
  end

let repo s1 s2 =
  Repo.[
    (* item ["protein_families"] s1#protein_families ; *)
    item ["assembly_metadata"] s1#assembly_metadata ;
    item ["assembly_list"] s1#assembly_list ;
  ]
