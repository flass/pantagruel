open Core
open Bistro
open Bistro_utils
open Tools

let dir_contents dir =
  Sys.readdir dir |> Array.to_list

let stage1 assembly_folder_path =
  let assembly_folder : _ dworkflow = Workflow.input assembly_folder_path in
  let assembly_list = assembly_list_of_folder assembly_folder in
  let assembly_stats = generate_assembly_stats assembly_folder in
  let assembly_info = assembly_info ~assembly_list in
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
    method assembly_info = assembly_info
    method assembly_metadata = assembly_metadata
    method protein_families = split_mmseqs_clustdb_fasta protein_families
    method all_unique_proteins = all_unique_proteins
  end


let stage2_pipeline fa =
  let alignment = Clustalo.clustalo fa in
  object
    method alignment = alignment
  end

let stage2 s1 ~assembly_folder family_folder =
  let orfan_fasta, family_fastas =
    match dir_contents family_folder with
    | h :: t -> h, List.take t 3 (* FIXME temp for test *)
    | [] -> assert false
  in
  let protein_alignment_folder =
    map_reduce_to_dir family_folder family_fastas ~f:Clustalo.clustalo
  in
  let cds_alignment_folder =
    extract_full_prot_and_cds_family_alignments
      ~protein_alignment_folder
      ~assembly_info:s1#assembly_info
      ~assembly_folder
      ~orfan_proteins:(Workflow.input (Filename.concat family_folder orfan_fasta))
  in
  object
    method cds_alignment_folder = cds_alignment_folder
    method protein_alignment_folder = protein_alignment_folder
  end

let repo s1 s2 =
  Repo.[
    item ["assembly_metadata"] s1#assembly_metadata ;
    item ["assembly_list"] s1#assembly_list ;
    item ["all_unique_proteins"] s1#all_unique_proteins ;
    item ["protein_families"] s1#protein_families ;
    item ["protein_alignments"] s2#protein_alignment_folder ;
    item ["cds_alignments"] s2#cds_alignment_folder ;
  ]
