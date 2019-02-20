open Core
open Bistro
open Bistro.Shell_dsl
open Bistro_bioinfo
open Bistro_utils

let _ = Workfl
(* configuration à fonctoriser FIXME *)
let email = "philippe.veber@univ-lyon1.fr"
let fam_prefix = "ABCDE"

let dir_contents dir =
  Sys.readdir dir |> Array.to_list

let tools_img = [ docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" () ]

let fetch_ncbi_taxonomy () =
  Workflow.shell ~descr:"fetch_ncbi_taxonomy" [
    cmd "bash" ~img:tools_img [
      file_dump (string Scripts.fetch_ncbi_taxonomy) ;
      string email ;
      dest ;
    ] ;
  ]

let assembly_list_of_folder (folder : _ dworkflow) =
  Workflow.shell ~descr:"assembly_list_of_folder" [
    cmd "ls" ~stdout:dest [
      string "-d" ;
      dep folder // "GC[AF]_*" ;
    ] ;
  ]

let collect_all_proteins (assembly_list : #text_file pworkflow)  =
  Workflow.shell ~descr:"collect_all_proteins" [
    cmd "bash" [
      file_dump (string Scripts.collect_all_proteins) ;
      dep assembly_list ;
      dest ;
    ] ;
  ]

let generate_assembly_stats assembly_folder =
  Workflow.shell ~descr:"generate_assembly_stats" [
    cmd "bash" ~img:tools_img [
      file_dump (string Scripts.generate_assembly_stats) ;
      dep assembly_folder ;
      dest ;
    ] ;
  ]

let extract_metadata_from_gbff ~assembly_list ~assembly_stats =
  Workflow.shell ~descr:"extract_metadata_from_gbff" [
    cmd "python" [
      file_dump (string Scripts.extract_metadata_from_gbff) ;
      opt "--assembly_folder_list" dep assembly_list ;
      opt "--add_assembly_info_dir" dep assembly_stats ;
      opt "--default_species_name" (string % quote ~using:'\'') "unclassified organism" ;
      opt "--output" ident dest ;
    ]
  ]

let dereplicate_fasta (fa : fasta pworkflow) : fasta pworkflow =
  Workflow.shell ~descr:"dereplicate_fasta" [
    cmd "python" ~img:tools_img [
      file_dump (string Scripts.dereplicate_fasta) ;
      dep fa ;
      dest ;
    ] ;
  ]

let assembly_info ~assembly_list =
  Workflow.shell ~descr:"allgenome_gff2db" [
    cmd "python" [
      file_dump (string Scripts.allgenome_gff2db) ;
      dep assembly_list ;
      dest ;
      dep (fetch_ncbi_taxonomy ()) // "scientific_names.dmp" ;
    ]
  ]

let split_mmseqs_clustdb_fasta fa =
  Workflow.shell ~descr:"split_mmseqs_clustdb_fasta" [
    cmd "python" ~img:tools_img [
      file_dump (string Scripts.split_mmseqs_clustdb_fasta) ;
      dep fa ;
      string (fam_prefix ^ "P") ;
      dest ;
    ]
  ]

let map_reduce_to_dir root files ~f =
  Workflow.shell ~descr:"map_reduce_to_dir" @@
  (
    mkdir_p dest ::
    List.map files ~f:(fun fn ->
        cmd "cp" [
          string "-p" ;
          dep (f (Workflow.input (Filename.concat root fn))) ;
          dest // fn ;
        ]
      )
  )

let extract_full_prot_and_cds_family_alignments
    ~protein_alignment_folder
    ~orfan_proteins
    ~assembly_folder
    ~assembly_info =
  Workflow.shell ~descr:"extract_full_prot_and_cds_family_alignments" [
    cmd "python" ~img:tools_img [
      file_dump (string Scripts.extract_full_prot_and_cds_family_alignments) ;
      dep protein_alignment_folder ;
      dep orfan_proteins ;
      dep assembly_info // "allproteins_info.tab" ;
      dep assembly_info // "allreplicons_info.tab" ;
      dep (Workflow.input assembly_folder) ;
      dest ;
      string fam_prefix ;
    ]
  ]

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
