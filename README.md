## Top Five Metal Coordination Geometry Scripts

This set of scripts runs the full set of analyses used in:

Yao S, Flight RM, Rouchka EC, Moseley HN. (2017) “Aberrant coordination geometries discovered in the most abundant metalloproteins.” Proteins: Structure, Function, and Bioinformatics, 85(5): 885-907 

## Citation and Licensing

Please see the individual packages for their license information, and the mainScripts/README.md and mainScripts/LICENSE for citing and license information of the metalloprotein analysis scripts.

## Installation

For a full script to full regenerate the data for the paper, including dependant packages, download the tarballed files at:

https://figshare.com/articles/Five_metal_manuscript/4229297

You can un-tar the downloaded file:

```
tar -xzpf fiveMetalManuscript.tar.gz
```

This creates the directory *fiveMetalManuscript*, that contains a set of directories:

  * atom2seq: R package that does alignment and functional enrichment
  * categoryCompare2: allows arbitrary enrichment calculations
  * mainScripts: the main scripts that do the analysis
  * output_metal_org: the original outputs generated for each and combined metal
  * metalPDB: supporting scripts in generating metalloprotein PDB list

## Prerequisites

To run it, you will need 32GB of Ram, and also need to install (on Linux):

```
perl::Math::MatrixReal
perl::Clone
perl::Thread::Queue
perl::JSON (optional for certain settings)
perl::Data::Dumper::Concise (optional for certain settings)

Interproscan (https://code.google.com/p/interproscan/)

libxml2-devel (for compiling some R packages from source)

R
R-devel
```

And then the following R packages:

```
Biostrings
Biobase
AnnotationDbi
Category
GO.db
annotate
clue
cluster
randomForest
ggplot2
gridExtra
```

In addition, a custom package also need to be installed, which are included. To install those, the easiest method is using `devtools`:

```
R
setwd("atom2seq")
devtools::install()

setwd("../categoryCompare2")
devtools::install()
```

## PDB download

You also need a local copy of the PDB.

`rsync -rlpt -v -z --delete --port=33444 \ > rsync.wwpdb.org::ftp_data/structures/divided/pdb/ ./pdb `

## Running analysis

### Counting metal sites

Use `./run.countMetal.bh` under the mainScripts directory to reproduce the full analysis.

### Main analysis

Finally, for actually running the scripts, you need to modify the alias for *interpro* and environment variables *output_dir* at the top of mainScripts/run.fullanalysis.bh.
Use `./run.fullanalysis.bh` under the mainScripts directory to reproduce the full analysis.

### 4-ligand Zn simulation
Use *./run.simulation.bh* under the mainScripts directory to reproduce the full analysis.

### Functional enrichmentment

Use `./run_enrichments.R` under the *mainScripts* directory to generate the enrichment results. The easiest way to run this is by:

```
cd mainScripts
R
source("run_enrichments.R")
```

The directories that are used may need to be changed.

### For manuscript figures

Run each R code in R console in the *supportingScripts* directory. Change input directory accordingly.

## Manuscript data

The original intermediate and final results generated for each and combined metal are in corresponding `output_*_org` directories.


Enjoy!
