options to create:

- step 0 (raw data downdload and grooming):  

- NCBI Assembly query:  
  - mandatory argument: `--taxid={int}`  
  - optional arguments for setting filters:  
```sh
# some categories of filters can be specified with option arguments to build the query
# all options are additive (build up the query string)
# these will be parsed as positive filters; they return a string like ' AND "filtername"[filter]'
--status={"latest"|"latest genbank"|"latest refseq"}
--assembly_level={"complete genome"|"chromosome"|"scaffold"|"contig"}
# these will be parsed as negative filters; they return a string like ' AND ("all"[filter] NOT "filtername"[filter])'
--exclude={"partial"|"anomalous"}
# can add raw filter specification
--filters='"filter1"[,"filter2"[, ...]]'
```
