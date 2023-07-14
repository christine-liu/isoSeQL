example_data contains 3 samples' files to be used as input to isoSeQL
Samples:
- NA18989
- NA19317
- NA19331

Each sample has a filtered classification file from running SQANTI3 (sample_filtered_classification.txt), a genePred file of all the isoforms from SQANTI3 (sample.genePred), a sample config file (sample_samp.config), and an experiment config file (sample_exp.config)


# Add to Database

In order to add these samples to a test database, run the following commands:
```
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA18989_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA18989.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA18989_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA18989_exp.config --db /path/to/output/example.db
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19317_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19317.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19317_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19317_exp.config --db /path/to/output/example.db
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19331_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19331.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19331_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19331_exp.config --db /path/to/output/example.db
```

# Query Database

Now that all three samples have been added to the database, they can be queried to make various plots/tables!

### What is currently in the database?
```
python /path/to/isoSeQL/isoSeQL_query.py expInfo --db /path/to/output/example.db --out /path/to/output/example_expInfo.txt
```
will generate /path/to/output/example_expInfo.txt, a text file including all sample/experiment information that has been added to the database.

For most queries, an "experiment list" is used as an input. This is a text file that lists the exp numbers to be queried, each on its own line. See examples /path/to/output/expList_all and /path/to/output/expList_1_3. Each of these lists can be used as queries to either include all 3 exps or only exps 1 and 3.

### Which isoforms are present in my samples?
```
python /path/to/isoSeQL/isoSeQL_query.py summary --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example
```
will generate two files:
- /path/to/output/example_commonJxn_allinfo_countMat.txt
- /path/to/output/example_ends_allinfo_countMat.txt
that keep track of isoform information (gene, category, subcategory, IEJ status) and corresponding counts for each exp. The difference between the two files is whether or not common junction isoforms (ends ignored) or isoforms with variable end coordinates are being examined.

### What proportion of isoforms/reads belong to the various structural categories (FSM, ISM, NIC, etc)?
```
python /path/to/isoSeQL/isoSeQL_query.py isoProp --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example
```
will generate 4 files:
- /path/to/output/example_isoPropPlot.pdf and /path/to/output/example_isoPropTable.txt, a stacked bar chart showing the proportion of ***isoforms*** that belong to each structural category and corresponding table
- /path/to/output/example_isoReadsPropPlot.pdf and /path/to/output/example_isoReadsPropPlot.txt, a stacked bar chart showing the proportion of ***reads*** supporting isoforms that belong to each structural category and corresponding table

The difference between these two plots is that one is not weighted by read count (isoform proportions) while the other is (read proportions)

### How can I see what these isoforms look like?
First make a names file to label the tracks in UCSC Genome Browser or IGV. This should be in the same order as the expList
```
python /path/to/isoSeQL/isoSeQL_query.py bed --db /path/to/output/example.db --exp /path/to/output/expList_1_3 --outPrefix /path/to/output/example --name nameList
```
will generate 1 file for every exp, /path/to/output/example_\[name\].bed. These files can be loaded into IGV or the UCSC Genome Browser. The isoforms are colored by structural category, and labeled by an identification number x_y_z, where x is the common isoform ID, y is the start coordinate, and z is the end coordinate.

### Which isoforms contain intra-exonic junctions (IEJs)?
IEJs can only be identified using [my altered version of SQANTI3](https://github.com/christine-liu/SQANTI3/tree/SQANTICL). This version of SQANTI3 provides additional annotations for NNC isoforms, identifying the different features that make them novel. IEJs were identified/defined in [this paper](https://www.nature.com/articles/s41586-018-0718-6) as junctions made up of a novel donor and novel acceptor that both occur within known annotated exons. These splice sites can both be within the same exon, in adjacent exons, or in non-adjacent distal exons.
```
python /path/to/isoSeQL/isoSeQL_query.py IEJtab --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example [--variable]
```
will generate /path/to/output/example_variableEnds_IEJs.txt if use --variable or /path/to/output/example_commonJxn_IEJs.txt otherwise. These tables contain the isoforms with IEJs and the corresponding read counts per exp. The difference between these files is whether or not you want to examine common junction isoforms or isoforms with variable ends (--variable).

### Can I generate input files (count matrix and gff) for tappAS?
```
python /path/to/isoSeQL/isoSeQL_query.py countMatrix --db /path/to/output/example.db --exp /path/to/output/expList_1_3 --outPrefix /path/to/output/example [--gene] [--variable]
```
there are 3 different file options, depending on which flags are used:
- no optional flags, /path/to/output/example_commonJxn_counts_matrix.txt, count matrix using common junction isoforms
- --gene, /path/to/output/example_gene_counts_matrix.txt, count matrix using ***genes*** (counts summed up over all isoforms of each gene)
- --variable, /path/to/output/example_variableEnds_counts_matrix.txt, count matrix using isoforms with variable ends

### Which isoforms are in common or unique to different exps?
```
python /path/to/isoSeQL/isoSeQL_query.py upset --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example [--top 20] [--variable]
```
The two optional flags allow you to switch from using common junction isoforms to isoforms with variable ends (--variable) and to change how many columns are displayed. The UpSet plot will be sorted from largest to smallest, and with the default settings, only the top 20 sets will be displayed. The larger the number, the longer it will take to graph.

output file options:
- --variable, /path/to/output/example_commonJxn_UpsetMatrix.txt and /path/to/output/example_varEnds_UpsetPlot.pdf, binary counts matrix used to create UpSet plot and resulting plot
- no flag, /path/to/output/example_varEnds_UpsetMatrix.txt and /path/to/output/example_commonJxn_UpsetPlot.pdf, binary counts matrix used to create UpSet plot and resulting plot

Bars are colored by structural categories of isoforms

### What does the relative expression of different known isoforms for a particular gene look like? (help display isoform switching)
First make a list of genes that you are interested in, each on their own line. /path/to/output/geneList
```
python /path/to/isoSeQL/isoSeQL_query.py FSM --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example --genes /path/to/output/geneList [--cutoff x]
```
if use --cutoff x, then the ***total*** number of reads for that gene needs to exceed x in order to plot. If not, an error message will inform you that "No samples exceed cutoff for \[gene\], please pick a different number" and no plot will be made.
if no cutoff is used, then data will be plotted regardless of read count.
For each gene in the geneList, a stacked bar plot indicating the proportion of reads from each known isoform for each exp will be generated. Numbers on the top of the bars will indicate total read count for that exp. /path/to/output/example_FSM_\[gene\].pdf

### What if none of these queries address what I'm looking for?
Lucky for you, your database can be loaded into python and you can write your own custom queries!
In python...
```
import sqlite3
import pandas as pd

conn=sqlite3.connect("/path/to/output/example.db")
c=conn.cursor()

df = pd.read_sql("INSERT QUERY HERE")
```
It's not always easy to come up with your own queries, so I recommend taking a look through the table set up to make sure you're getting the numbers that you expect