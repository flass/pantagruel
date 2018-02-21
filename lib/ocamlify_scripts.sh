ocamlify $(for f in `ls scripts/*[^~]`; do echo --var-string $(basename ${f%.*}) $f; done | xargs) --output scripts.ml
