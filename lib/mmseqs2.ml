open Bistro.EDSL

let env = docker_image ~account:"pveber" ~name:"mmseqs2" ~tag:"c7a89" ()

let subcmd subcmd args =
  cmd "mmseqs" ~env (string subcmd :: args)

let createdb fa =
  workflow ~descr:"mmseqs2.createdb" [
    mkdir_p dest ;
    subcmd "createdb" [
      dep fa ;
      dest // "mmseqs-db" ;
    ];
  ]

let cluster db =
  workflow ~descr:"mmseq2.cluster" [
    mkdir_p dest ;
    subcmd "cluster" [
      dep db // "mmseqs-db" ;
      dest // "mmseqs-clustering" ;
      tmp ;
    ]
  ]

let createseqfiledb db clusters =
  workflow ~descr:"mmseq2.createseqfiledb" [
    subcmd "createseqfiledb" [
      dep db // "mmseqs-db" ;
      dep clusters // "mmseqs-clustering" ;
      dest ;
    ] ;
  ]
