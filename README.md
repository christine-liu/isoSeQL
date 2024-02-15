# isoSeQL
Database tool for comparing across many Iso-Seq runs analyzed through SQANTI3

SQANTI3: https://github.com/ConesaLab/SQANTI3

## Latest Updates
(July, 13, 2023) isoSeQL v1.0.0 is now live!


## How to cite isoSeQL
The isoSeQL manuscript is currently in preparation to be submitted. A bioRxiv version may be available soon

## Installation
isoSeQL was developed using python 3.9.6

Dependencies:

R (v 4.1.1)

R libraries optparse and ComplexUpset

Python libraries numpy, pandas, plotly

Using anaconda3
```
conda create --name isoSeQL python=3.9.6
conda activate isoSeQL
conda install -c conda-forge r-base=4.1.1
conda install -c plotly plotly
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c conda-forge python-kaleido
conda install -c conda-forge r-devtools
```
In R, install optparse and ComplexUpset
```
install.packages("optparse")
devtools::install_github("krassowski/complex-upset")
```
I've included a yml file to help build the conda environment as well
```
conda env create -f isoSeQL.yml
```

## Terminology
_common junction isoforms_ - structure information stored in the **isoform** table, consolidates all isoforms with the same junctions, regardless of start and end coordinate, counts stored in **counts** table

_variable ends isoforms/isoforms with variable ends_ - structure information linked to common junction isoform in **isoform** table, start/end coordinates stored in **isoform_ends** table, counts stored in **ends_counts** table

_experiment (exp)_ vs _sample_ - I wanted to be able to link together two different experiments (run at different times or together) that used the same sample material. **sampleData** stores all the information about the cell line, tissue, etc being used as source material. **exp** tracks the software versions, reference versions, and when the experiment was performed. Currently queries can be performed to examine isoforms by _exp_ ID (integers). Future implementations will include the ability to group exp IDs together if they come from the same sample or same sample group.

## How to use isoSeQL
See example/example.md or the Wiki for example commands that walk you through adding samples to a database and all the current built-in queries (as of January 2024) !

## On the way
Additional queries will be added as I write new functions!
