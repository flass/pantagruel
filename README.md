![repas]

# required softwares
- **MMseqs2/Linclust** for homologous sequence clustering  
  (Install from [source code](https://github.com/soedinglab/MMseqs2); last tested version https://github.com/soedinglab/MMseqs2/commit/c92411b91175a2362554849b8889a5770a1ae537)

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
  (Install from [source code](https://github.com/ssolo/ALE); version used and recommended: 0.4; notably depends on [Bio++ libs](https://github.com/BioPP) (v2.2.0))
  

- R (version 3, >=3.2.3 recommended) + packages:
  - ape
  - phytools
  - DBI, RSQLite
  
- Python (version 2.7, >=2.7.13 recommended) + packages:
  - [sqlite3](https://docs.python.org/2/library/sqlite3.html) (standard package in Python 2.7)
  - [scipy/numpy](https://www.scipy.org/scipylib/download.html)
  - [tree2](https://github.com/flass/tree2)
  - [BioPython](http://biopython.org/wiki/Download)

- Other software:
  - [sqlite3](https://www.sqlite.org) (available as a Debian package)
  - [LFTP](https://lftp.yar.ru/get.html) (available as a Debian package)

-------------

Two version of the pipeline are distributed:  

- a script version, which source code is adaptable and can be deployed on high-performance computing (HPC) "cluster" Linux systems;  

- a pre-compiled Docker image that can be deployed on pretty much any platform, including swarms of virtual machines (VMs). The latter version was implemented using P. Veber's [Bistro](https://github.com/pveber/bistro) framework.


[repas]: https://github.com/flass/pantagruel/blob/master/Pantagruels_childhood.jpg
