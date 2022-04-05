#!/usr/bin/env python


#07/29/2021
import sqlite3
import argparse


#Generate database and its tables
conn=sqlite3.connect("testDatabase.db") #makes database if doesn't already exist
c=conn.cursor()


c.execute("CREATE TABLE isoform (id INTEGER PRIMARY KEY, chr TEXT, strand TEXT, junctions TEXT, gene TEXT, iso_exons INTEGER, subcategory TEXT, canonical TEXT, IEJ INTEGER, category TEXT, UNIQUE(chr, strand, junctions))")
c.execute("CREATE TABLE isoform_ends (id INTEGER PRIMARY KEY, isoform_id INTEGER, chr TEXT, start INTEGER, end INTEGER, exp INTEGER,read_count INTEGER, FOREIGN KEY(isoform_id) REFERENCES isoform(id), FOREIGN KEY(exp) REFERENCES exp(id)), UNIQUE(isoform_id, exp, start, end)")
c.execute("CREATE TABLE counts (id INTEGER PRIMARY KEY, isoform_id INTEGER, exp INTEGER, read_count INTEGER, FOREIGN KEY(isoform_id) REFERENCES isoform(id),FOREIGN KEY(exp) REFERENCES exp(id))")
c.execute("CREATE TABLE exp (id INTEGER PRIMARY KEY, patient_id INTEGER, seq_date DATE, vMap TEXT, vReference TEXT, vAnnot TEXT, vLima TEXT, vCCS TEXT, vIsoseq3 TEXT, vCupcake TEXT, vSQANTI TEXT, FOREIGN KEY(patient_id) REFERENCES sampleData(id))")
c.execute("CREATE TABLE sampleData (id INTEGER PRIMARY KEY, patient TEXT, tissue TEXT, disease TEXT, RIN REAL, age INTEGER, UNIQUE(patient, tissue, disease, age))")
c.execute("CREATE TABLE PBID (id INTEGER PRIMARY KEY, PBID TEXT, exp INTEGER, isoform_id INTEGER, FOREIGN KEY(exp) REFERENCES exp(id), FOREIGN KEY(isoform_id) REFERENCES isoform(id)), UNIQUE(PBID, exp, isoform_id)")


#Check if sampleData added - has this sample been analyzed with iso-seq before?
#parse information about sample from sampleData file?
#check if sample already exists, if not then add and get sampleID, if so then get sampleID
c.execute('SELECT id FROM sampleData WHERE patient = ? AND tissue = ? and disease = ? and age = ?', (patient, tissue, disease, age,))
sampleID = c.fetchall()
if len(sampleID) == 0:
	c.execute('INSERT INTO sampleData(patient, tissue, disease, RIN, age) VALUES (?,?,?,?,?)', (patient, tissue, disease, RIN, age,))
	sampleID = c.lastrowid
else:
	sampleID=sampleID[0][0]

#Add exp info (each sample gets analyzed individually, so will refer to as an experiment b/c can run the same sample multiple times in different runs)
#parse info about exp from an exp data file
c.execute("INSERT INTO exp(patient_id, seq_date, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (sampleID, date, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI,))
expID = c.lastrowid

##insert parsing of genepred and classification files into dictionaries so that we can just refer to those
classif = parse_classification(ARGS)
genePred = parse_genePred(ARGS)

observedIsoIDs = []
for iso in classif.keys():
	c.execute('SELECT id FROM isoform WHERE chr = ? AND strand = ? AND junctions = ?', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions)))
	isoID = c.fetchall() ##check if this is in a weird format that need to reconfigure
	if len(isoID) == 0:
		c.execute('INSERT INTO isoform(chr, strand, junctions, gene, iso_exons, subcategory, canonical, IEJ, category) VALUES (?,?,?,?,?,?,?,?,?)', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions), classif[iso].gene, classif[iso].tx_exons, classif[iso].subcat, classif[iso].canonical, classif[iso].IEJ, classif[iso].cat,))
		isoID=c.lastrowid
	else:
		isoID=isoID[0][0]
	observedIsoIDs.append(isoID)
	c.execute('INSERT INTO isoform_ends(isoform_id, chr, start, end, exp, read_count) VALUES (?,?,?,?,?,?)', (isoID, genePred[iso].chrom, genePred[iso].start, genePred[iso].end, expID, classif[iso].count))
	c.execute('INSERT INTO PBID(PBID, exp, isoform_id) VALUES (?,?,?)', (iso, expID, isoID))
#sum up counts by querying each isoform ID that was observed and adding to counts table
for id in observedIsoIDs:
	c.execute('SELECT SUM(read_count) FROM isoform_ends WHERE isoform_id = ? AND read_count != "NA"', (id,))
	sum_counts = c.fetchall()
	c.execute('INSERT INTO counts(isoform_id, exp, read_count) VALUES (?,?,?)', (id, expID, sum_counts))


conn.commit()





