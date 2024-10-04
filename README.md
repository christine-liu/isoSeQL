![isoSeQL logo](images/isoSeQL_logo.png)
# isoSeQL
Database tool for comparing across many Iso-Seq runs analyzed through SQANTI3

SQANTI3: https://github.com/ConesaLab/SQANTI3

# Latest Updates
(October 7, 2024) isoSeQL v1.0.0 official release tagged on Github

(July 13, 2023) isoSeQL v1.0.0 is now live!

# Motivation (Origin Story)
isoSeQL was created to work with the output of SQANTI2 (now SQANTI3). At the time that I started working with long-read Iso-Seq data (\~end of 2019), I followed the recommended workflow for analyzing PacBio Iso-Seq data that ended in annotating the high quality isoforms with SQANTI2. I was working with multiple samples that I wanted to compare to each other, but quickly realized that the isoform ID #s (PB.X.Y) in each sample couldn't be matched. Since I couldn't find a tool to do this for me, and I really needed to figure out how to analyze the data, I solved this issue by labeling all the HiFi reads from each sample, concatenating the read files together, running the analysis pipeline (isoseq > minimap2 > cDNA_cupcake > SQANTI2) on that single "mega sample", and deconvoluting which reads from which samples supported each isoform after everything was annotated. This method 100% works, but does not scale well with increased sample numbers or especially as the throughput of SMRTcells has increased over the years. I wanted to figure out a way to work directly with the SQANTI2/3 output files to consolidate and compare across multiple samples.

I want to acknowledge that many other tools are being developed, and the long-read isoform sequencing field is constantly progressing. There are tools that work with data upstream of SQANTI3 to create transcript models, annotate, and quantify isoforms, etc. Not only did many of these tools not exist when I was first figuring out how to work with Iso-Seq data from multiple samples, but I had also made some changes to SQANTI3 itself in order to identify novel isoform features that I wanted to examine more closely. So I needed to find a way to make comparisons using the SQANTI3 outputs specifically instead of trying out a new tool/method. This is a pretty specific usage instance for isoSeQL, but I've seen many people experience the same confusion I did upon realizing that they've created SQANTI3-annotated isoform files than cannot be compared, so I hope that this tool can be useful to many others!

Happy analyzing!

# How to cite isoSeQL
The isoSeQL manuscript is currently in preparation to be submitted. A bioRxiv version may be available soon

# Installation
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

# Terminology
<!-- _common junction isoforms_ - structure information stored in the **isoform** table, consolidates all isoforms with the same junctions, regardless of start and end coordinate, counts stored in **counts** table

_variable ends isoforms/isoforms with variable ends_ - structure information linked to common junction isoform in **isoform** table, start/end coordinates stored in **isoform_ends** table, counts stored in **ends_counts** table

 -->_experiment (exp)_ vs _sample_ - I wanted to be able to link together two different experiments (run at different times or together) that used the same sample material. **sampleData** stores all the information about the cell line, tissue, etc being used as source material. **exp** tracks the software versions, reference versions, and when the experiment was performed. Currently queries can be performed to examine isoforms by _exp_ ID (integers). Future implementations will include the ability to group exp IDs together if they come from the same sample or same sample group.

# How to use isoSeQL

## Upstream bulk Iso-Seq data analysis (before input into isoSeQL)
I've listed out the steps required to process samples before they are analyzed using isoSeQL. This workflow generally follows the [recommendations from PacBio](https://isoseq.how/) as of September 2024 and requires additional packages to be installed. Please refer to their workflow for additional details.

### Generating HiFi reads
If you're starting with sample_subreads.bam, then you'll want to use the [ccs tool](https://ccs.how) to generate HiFi reads. Simply put, this step is generating a consensus read from the multiple polymerase passes across the cDNA. 
```
ccs sample_subreads.bam --report-file /outdir/sample_isoSeqccs_report.txt /outdir/sample_isoSeqccs.bam

```
This outputs sample_isoSeqccs.bam (bam of HiFi reads) and sample_isoSeqccs_report.txt which summarizes how many ZMWs passed or failed.

### Removing Iso-Seq primers and orienting sequences
[lima](https://lima.how/) is used to identify and remove primer sequences using a specialized --isoseq mode that recognizes outputs reads specifically with the desired/expected asymmetric primer combination (5' and 3'). Lima will also orient sequences in the 5' to 3' direction. 
```
lima /outdir/sample_isoSeqccs.bam /path/to/primers.fa /outdir/sample_isoSeqccs_demux.bam --isoseq --peek-guess 

```
This outputs sample_isoSeqccs_demux.bam (bam of full-length, properly oriented reads)

### Refine full-length reads
[isoseq](https://isoseq.how) refine is used to trim polyA tails and remove concatemeric sequences
```
isoseq refine --require-polya /outdir/sample_isoSeqccs_demux.Primer_5p--Primer_3p.bam /path/to/primers.fa /outdir/sample_isoSeqccs_demux_flnc.bam
```
This outputs sample_isoSeqccs_demux_flnc.bam (bam of full-length, non-concatemeric (FLNC) reads)

### Cluster similar sequences
[isoseq](https://isoseq.how) cluster2 generates clusters of isoforms with the support of at least 2 FLNC reads
```
isoseq cluster2 /outdir/sample_isoSeqccs_demux_flnc.bam /outdir/sample_isoSeqccs_demux_flnc_cluster.bam
```
This outputs sample_isoSeqccs_demux_flnc_cluster.bam, sample_isoSeqccs_demux_flnc_cluster.hq.fasta.gz (isoforms with predicted accuracy >= 0.99), sample_isoSeqccs_demux_flnc_cluster.lq.fasta.gz (isoforms with predicted accuracy < 0.99)

### Map to the reference genome
[pbmm2](https://github.com/PacificBiosciences/pbmm2) is a wrapper for minimap2 that has sets of recommended parameters that support PacBio data alignment.
```
pbmm2 align --preset ISOSEQ -G 2500000 --sort /path/to/reference.fa /outdir/sample_isoSeqccs_demux_flnc_cluster.bam /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped.bam
```
This outputs sample_isoSeqccss_demux_flnc_cluster_mapped.bam

### Collapse redundant transcripts
[isoseq](https://isoseq.how) collapse is used to generate a set of unique isoforms. Previously (cDNA_Cupcake)[https://github.com/Magdoll/cDNA_Cupcake] was used for this step, and while it is deprecated, it still works well, but significantly slower than isoseq collapse.
```
isoseq collapse /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped.bam /outdir/sample_isoSeqccs/demux_flnc.bam /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed.gff
```
This outputs sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed.gff and sample_isoSeq_demux_flnc_cluster_mapped_collapsed.flnc_count.txt.

### Annotate and filter isoforms using SQANTI3
As mentioned in the Motivation/Origin Story, I made some modifications to SQANTI3 to add some additional annotation details about novel not in catalog (NNC) isoforms - specifically which features present in the isoform make it deviate from the known, reference isoforms. These modifications and features are described in a paper that is currently under review, which I will link upon publication. This modified version of SQANTI3 is available [here](https://github.com/christine-liu/SQANTI3/tree/SQANTICL). 
```
python /path/to/SQANTI3/sqanti3_qc.py /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed.gff /path/to/annotation.gtf /path/to/reference.fa --genename --report skip --fl_count /outdir/sample_isoSeq_demux_flnc_cluster_mapped_collapsed.flnc_count.txt

```
This outputs numerous files which can then be used to filter for high-confidence isoforms (also using SQANTI3). Due to the modifications that I made, I continue to use the rules filter. I skip the report generation so that that part doesn't error out due to the modifications that I made. Supplying the \*flnc_count.txt file populates a column in the \*classification.txt file with the number of reads supporting that isoform.
```
python /path/to/SQANTI3/sqanti3_RulesFilter.py /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed_classification.txt /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed_corrected.fasta /outdir/sample_isoSeqccs_demux_flnc_cluster_mapped_collapsed.gff --report skip
```
This outputs filtered versions of several files, including the classfication file that will be used as an input to isoSeQL.

## Analysis with isoSeQL
The [example page](example/example.md) and the Wiki also have example commands that walk you through adding samples to a database and all the current built-in queries (as of September 2024)!


# On the way
Additional queries will be added as I write new functions!
