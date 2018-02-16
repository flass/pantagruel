![repas]

# required softwares
- **MMseqs2** for homologous sequence clustering  
  (Install from [source code](https://github.com/soedinglab/MMseqs2)

- **Clustal Omega** for homologous sequence alignment  
  (Install from [source code](http://www.clustal.org/omega/) or *clustalo* debian package; version used and recommended: 1.2.1)  
  - \[ future development: consider using [FAMSA](http://sun.aei.polsl.pl/REFRESH/famsa) \]

- **PAL2NAL** for reverse tanslation of protein sequence alignments into CDS alignments  
  ([Perl source code](http://www.bork.embl.de/pal2nal/))

- **RAxML** for species tree and initial (full) gene tree estimation  
  (Install from [source code](https://github.com/stamatak/standard-RAxML) or *raxml* debian package; version used and recommended: 8.2.9)  
  - \[ future development: consider using RAxML-NG (Install from [source code](https://github.com/amkozlov/raxml-ng)) \]

- **MrBayes** for secondary estimation of (collapsed) gene trees  
  (Install from [source code](http://mrbayes.sourceforge.net/) or *mrbayes* and *mrbayes-mpi* debian packages; version used and recommended: 3.2.6)  
  - \[ future development: consider using [RevBayes](http://revbayes.github.io/) \]

- **MAD** for species tree rooting  
  ([R source code](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad-r-tar.gz))

- **ALE/xODT** for gene tree / species tree reconciliation  
  (Install from [source code](https://github.com/ssolo/ALE))
  

- R packages:
  - ape
  - phytools
  
- Python packages:
  - [psycopg2](https://pypi.python.org/pypi/psycopg2)
  - [scipy/numpy](https://www.scipy.org/scipylib/download.html)
  - [tree2](https://github.com/flass/tree2)
  - [BioPython](http://biopython.org/wiki/Download)


[repas]: https://github.com/flass/pantagruel/blob/master/Pantagruels_childhood.jpg
