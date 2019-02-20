open Core
open Bistro
open Bistro.Shell_dsl

let img = [ docker_image ~account:"pveber" ~name:"clustalo" ~tag:"1.2.4" () ]

let clustalo fa =
  Workflow.shell ~descr:"clustalo" [
    cmd "clustalo" ~img [
      opt "-i" dep fa ;
      opt "-o" ident dest ;
    ] ;
  ]
