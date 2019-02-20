open Core
open Bistro
open Bistro.Shell_dsl
open Bistro_bioinfo

(* configuration à fonctoriser FIXME *)
let email = "philippe.veber@univ-lyon1.fr"
let fam_prefix = "ABCDE"

let img = [ docker_image ~account:"pveber" ~name:"pantagruel-tools" ~tag:"latest" () ]

let fetch_ncbi_taxonomy () =
  Workflow.shell ~descr:"fetch_ncbi_taxonomy" [
    cmd "bash" ~img [
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
    cmd "bash" ~img [
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
    cmd "python" ~img [
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
    cmd "python" ~img [
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
    cmd "python" ~img [
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
