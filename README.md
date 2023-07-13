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

R libraries optparse and UpSetR

Python libraries numpy, pandas, plotly

## Terminology
The isoSeQL SQLite database is made up of several different tables.

In order to unify isoforms IDs across different samples and deal with the challenge of isoforms having variable start and end coordinates due to fragmentation or other factors, isoform structural information is split between two tables **isoform** and **isoform_ends**. **isoform** keeps track of the chromosome, strand, junctions (encoded as a string '(a,b),(c,d),(e,f)'), gene, number of exons, structural category (FSM,ISM,NIC,etc), subcategory (if applicable), canonical (T/F for use of canonical splice sites), IEJ (T/F for containing 1+ IEJ). These isoforms are referred to as "common junction isoforms" b/c it collapses together isoforms that have the exact same junctions but different start and/or end coordinates and treats them as a single isoform. **isoform_ends** is linked to the corresponding common junction isoform in **isoform** but additionally keeps track of the start and end coordinates, exon sizes, and exon starts. These isoforms are referred to as "isoforms with variable ends". From these two tables alone, it's possible to generate a file (like a gff) with all the structural information for each isoform.

Two tables **exp** and **sampleData** keep track of metadata associated with each addition to the database. The logic behind keeping track of experiments (exp) separate from samples is that certain tissue samples, cells, etc may be used for multiple runs. For example, if I use a piece of brain tissue for an experiment and then re-use it a year later, I want to be able to link those two experiments together for comparison b/c I used the same tissue sample. I still want to be able to treat the two sequencing runs separately though, and if they were run at different times, it's important to keep track of any software or reference versions that were used for the analysis of each run. The **sampleData** table keeps track of the sample name, tissue, disease status, age, and sex. Its linked to various entries in the **exp** table which additionally keeps track of RIN (RNA Integrity Number, measures quality/fragmented-ness of RNA sample, correlates with sequencing quality metrics), date of experiment, platform (ONT, PacBio Sequel I/II/IIe/Revio), method (bulk, single-cell), mapping software version, reference genome, reference annotation, lima version (demux software), ccs version (circular consensus sequencing software), isoseq3 version, cupcake version (previously used for collapsing isoforms, but now isoseq3 can be used for collapsing), SQANTI version, experiment name. Currently queries can be performed to examine isoforms by _exp_ ID (integers). Future implementations will include the ability to group exp IDs together if they come from the same sample or same sample group.


## How to use isoSeQL
See example/example.md for example commands that walk you through adding samples to a database and all the current built-in queries (as of July 2023) !

## On the way
isoSeQL has been adapted to work with single-cell long-read data (not part of publication). It needs optimization, and a tutorial will follow when that is completed.
