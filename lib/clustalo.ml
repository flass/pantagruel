open Core
open Bistro.EDSL

let env = docker_image ~account:"pveber" ~name:"clustalo" ~tag:"1.2.4" ()

let clustalo fa =
  workflow ~descr:"clustalo" [
    cmd "clustalo" ~env [
      opt "-i" dep fa ;
      opt "-o" ident dest ;
    ] ;
  ]
