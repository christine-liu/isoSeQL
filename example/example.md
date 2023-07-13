example_data contains 3 samples' files to be used as input to isoSeQL
Samples:
- NA18989
- NA19317
- NA19331

Each sample has a filtered classification file from running SQANTI3 (sample_filtered_classification.txt), a genePred file of all the isoforms from SQANTI3 (sample.genePred), a sample config file (sample_samp.config), and an experiment config file (sample_exp.config)


#Add to Database

In order to add these samples to a test database, run the following commands:

```
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA18989_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA18989.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA18989_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA18989_exp.config --db /path/to/output/example.db
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19317_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19317.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19317_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19317_exp.config --db /path/to/output/example.db
python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19331_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19331.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19331_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19331_exp.config --db /path/to/output/example.db
```

#Query Database

Now that all three samples have been added to the database, they can be queried to make various plots/tables!

###What is currently in the database?
```
python /path/to/isoSeQL/isoSeQL_query.py expInfo --db /path/to/output/example.db --out /path/to/output/example_expInfo.txt
```
will generate /path/to/output/example_expInfo.txt, a text file including all sample/experiment information that has been added to the database.

For most queries, an "experiment list" is used as an input. This is a text file that lists the exp numbers to be queried, each on its own line. See examples /path/to/output/expList_all and /path/to/output/expList_1_3. Each of these lists can be used as queries to either include all 3 exps or only exps 1 and 3.

###What proportion of isoforms/reads belong to the various structural categories (FSM, ISM, NIC, etc)?

```
python /path/to/isoSeQL/isoSeQL_query.py isoProp --db /path/to/output/example.db --exp /path/to/output/expList_all --outPrefix /path/to/output/example
```
will generate 4 files:
- /path/to/output/example_isoPropPlot.pdf and /path/to/output/example_isoPropTable.txt, a stacked bar chart showing the proportion of ***isoforms*** that belong to each structural category and corresponding table
- /path/to/output/example_isoReadsPropPlot.pdf and /path/to/output/example_isoReadsPropPlot.txt, a stacked bar chart showing the proportion of ***reads*** that belong to each structural category and corresponding table

The difference between these two plots is that one is not weighted by read count (isoform proportions) while the other is (read proportions)

