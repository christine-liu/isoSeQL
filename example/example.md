example_data contains 3 samples' files to be used as input to isoSeQL
Samples:
- NA18989
- NA19317
- NA19331

Each sample has a filtered classification file from running SQANTI3 (sample_filtered_classification.txt), a genePred file of all the isoforms from SQANTI3 (sample.genePred), a sample config file (sample_samp.config), and an experiment config file (sample_exp.config)


===Add to Database===
In order to add these samples to a test database, run the following commands:

python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA18989_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA18989.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA18989_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA18989_exp.config --db /path/to/output/example.db

python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19317_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19317.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19317_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19317_exp.config --db /path/to/output/example.db

python /path/to/isoSeQL/isoSeQL_run.py --classif /path/to/isoSeQL/example/example_data/NA19331_filtered_classification.txt --genePred /path/to/isoSeQL/example/example_data/NA19331.genePred --sampleConfig /path/to/isoSeQL/example/example_data/NA19331_samp.config --expConfig /path/to/isoSeQL/example/example_data/NA19331_exp.config --db /path/to/output/example.db

===Query Database===
Now that all three samples have been added to the database, they can be queried to make various plots/tables!


