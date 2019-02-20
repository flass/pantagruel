open Bistro
open Bistro.Shell_dsl

let img = [ docker_image ~account:"pveber" ~name:"mmseqs2" ~tag:"c7a89" () ]

let subcmd subcmd args =
  cmd "mmseqs" ~img (string subcmd :: args)

let createdb fa =
  Workflow.shell ~descr:"mmseqs2.createdb" [
    mkdir_p dest ;
    subcmd "createdb" [
      dep fa ;
      dest // "mmseqs-db" ;
    ];
  ]

let cluster db =
  Workflow.shell ~descr:"mmseq2.cluster" [
    mkdir_p dest ;
    subcmd "cluster" [
      dep db // "mmseqs-db" ;
      dest // "mmseqs-clustering" ;
      tmp ;
    ]
  ]

let createseqfiledb db clusters =
  Workflow.shell ~descr:"mmseq2.createseqfiledb" [
    subcmd "createseqfiledb" [
      dep db // "mmseqs-db" ;
      dep clusters // "mmseqs-clustering" ;
      dest ;
    ] ;
  ]
