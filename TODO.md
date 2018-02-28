find a story behind the name! or a derivation of the name.  

- PANTEGRUALE for "PANgenome inTEGRation Using ALE" ?

-------------

options to create:

- step 0 (raw data downdload and grooming):  

- NCBI Assembly query:  
  - mandatory argument: `--taxid={int}`  
  - optional arguments for setting filters:  
```bash
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

- possibility to build image with different instruction sets
  - need to detect supported instruction sets: 
```bash
# find newest arch support
for inset in avx2 avx sse4_2 sse4_1 sse3 sse2 ; do 
  if [ ! -z "$(grep $inset /proc/cpuinfo)" ] ; then
    echo $inset
    break
  fi
done
# build tools (MMSeqs2, RAxML, ...) with the appropriate instruction set

```
