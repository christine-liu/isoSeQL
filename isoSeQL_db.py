#!/usr/bin/env python

import sqlite3
import csv

def make_db(database):
	conn=sqlite3.connect(database)
	c=conn.cursor()
	c.execute("CREATE TABLE isoform (id INTEGER PRIMARY KEY, chr TEXT, strand TEXT, junctions TEXT, gene TEXT, iso_exons INTEGER, subcategory TEXT, canonical TEXT, IEJ INTEGER, category TEXT, UNIQUE(chr, strand, junctions, gene))")
	c.execute("CREATE TABLE isoform_ends (id INTEGER PRIMARY KEY, isoform_id INTEGER, chr TEXT, start INTEGER, end INTEGER, exp INTEGER, read_count INTEGER, ex_sizes TEXT, ex_starts TEXT, FOREIGN KEY(isoform_id) REFERENCES isoform(id), FOREIGN KEY(exp) REFERENCES exp(id), UNIQUE(isoform_id, exp, start, end))")
	c.execute("CREATE TABLE counts (id INTEGER PRIMARY KEY, isoform_id INTEGER, exp INTEGER, read_count INTEGER, FOREIGN KEY(isoform_id) REFERENCES isoform(id),FOREIGN KEY(exp) REFERENCES exp(id), UNIQUE(isoform_id, exp))")
	c.execute("CREATE TABLE exp (id INTEGER PRIMARY KEY, sample_id INTEGER, RIN REAL, seq_date DATE, platform TEXT, vMap TEXT, vReference TEXT, vAnnot TEXT, vLima TEXT, vCCS TEXT, vIsoseq3 TEXT, vCupcake TEXT, vSQANTI TEXT, exp_name TEXT, FOREIGN KEY(sample_id) REFERENCES sampleData(id) UNIQUE(sample_id, RIN, seq_date, platform, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI, exp_name))")
	c.execute("CREATE TABLE sampleData (id INTEGER PRIMARY KEY, sample_name TEXT, tissue TEXT, disease TEXT, age INTEGER, sex TEXT,UNIQUE(sample_name, tissue, disease, age, sex))")
	c.execute("CREATE TABLE PBID (id INTEGER PRIMARY KEY, PBID TEXT, exp INTEGER, isoform_id INTEGER, FOREIGN KEY(exp) REFERENCES exp(id), FOREIGN KEY(isoform_id) REFERENCES isoform(id), UNIQUE(PBID, exp, isoform_id))")
	c.execute("CREATE TABLE txID (id INTEGER PRIMARY KEY, tx TEXT, exp INTEGER, isoform_id INTEGER, gene TEXT, FOREIGN KEY(exp) REFERENCES exp(id), FOREIGN KEY(isoform_id) REFERENCES isoform(id), UNIQUE(tx, exp, isoform_id))")
	conn.commit()
	conn.close()

def addSampleData(database, sampleConfig): #returns sampleID for use later
	conn=sqlite3.connect(database)
	c=conn.cursor()

	file=open(sampleConfig, "rt")
	sampleRead=csv.DictReader(file)
	for row in sampleRead:
		c.execute('SELECT id FROM sampleData WHERE sample_name = ? AND tissue = ? and disease = ? and age = ? and sex = ?', (row['sample_name'], row['tissue'], row['disease'], row['age'], row['sex']))
		sampleID = c.fetchall()
		if len(sampleID) == 0:
			c.execute('INSERT INTO sampleData(sample_name, tissue, disease, age, sex) VALUES (?,?,?,?,?)', (row['sample_name'], row['tissue'], row['disease'], row['age'], row['sex'],))
			sampleID = c.lastrowid
		else:
			sampleID=sampleID[0][0]
	conn.commit()
	conn.close()
	return sampleID

def addExpData(database, expConfig, sampleID): #returns expID for use later
	conn=sqlite3.connect(database)
	c=conn.cursor()

	file=open(expConfig, "rt")
	expRead=csv.DictReader(file)
	for row in expRead:
		c.execute('SELECT id FROM exp WHERE sample_id = ? AND RIN = ? AND seq_date = ? AND platform = ? AND vMap = ? AND vReference = ? AND vAnnot = ? AND vLima = ? AND vCCS = ? AND vIsoseq3 = ? AND vCupcake = ? AND vSQANTI = ? AND exp_name = ?', (sampleID, row['RIN'], row['date'], row['platform'], row['vMap'], row['vReference'], row['vAnnot'], row['vLima'], row['vCCS'], row['vIsoseq3'], row['vCupcake'], row['vSQANTI'], row['exp_name'],))
		expID = c.fetchall()
		if len(expID) == 0:
			c.execute("INSERT INTO exp(sample_id, RIN, seq_date, platform, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI, exp_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (sampleID, row['RIN'], row['date'], row['platform'], row['vMap'], row['vReference'], row['vAnnot'], row['vLima'], row['vCCS'], row['vIsoseq3'], row['vCupcake'], row['vSQANTI'], row['exp_name'],))
			expID=c.lastrowid
		else:
			expID = "Already in database"
	conn.commit()
	conn.close()
	return expID

def addIsoforms(database, classif, genePred, expID):
	conn=sqlite3.connect(database)
	c=conn.cursor()
	observedIsoIDs = set()
	for iso in classif.keys():
		c.execute('SELECT id FROM isoform WHERE chr = ? AND strand = ? AND junctions = ? and gene = ?', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions), classif[iso].gene))
		isoID = c.fetchall()
		if len(isoID) == 0:
			c.execute('INSERT INTO isoform(chr, strand, junctions, gene, iso_exons, subcategory, canonical, IEJ, category) VALUES (?,?,?,?,?,?,?,?,?)', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions), classif[iso].gene, classif[iso].tx_exons, classif[iso].subcat, classif[iso].canonical, classif[iso].IEJ, classif[iso].cat,))
			isoID=c.lastrowid
		else:
			isoID=isoID[0][0]
		observedIsoIDs.add(isoID)
		c.execute('INSERT INTO isoform_ends(isoform_id, chr, start, end, exp, read_count, ex_sizes, ex_starts) VALUES (?,?,?,?,?,?,?,?)', (isoID, genePred[iso].chrom, genePred[iso].start, genePred[iso].end, expID, classif[iso].count, genePred[iso].exSizes, genePred[iso].exBedStarts))
		c.execute('INSERT INTO PBID(PBID, exp, isoform_id) VALUES (?,?,?)', (iso, expID, isoID))
		if classif[iso].cat == "full-splice_match" or classif[iso].cat=="incomplete-splice_match":
			c.execute('INSERT OR IGNORE INTO txID(isoform_id, exp, tx, gene) VALUES (?,?,?,?)', (isoID, expID, classif[iso].transcript, classif[iso].gene))
	for id in observedIsoIDs:
		c.execute('SELECT SUM(read_count) FROM isoform_ends WHERE isoform_id = ? AND read_count != "NA" AND exp = ?', (id,expID))
		sum_counts = int(c.fetchall()[0][0])
		c.execute('INSERT INTO counts(isoform_id, exp, read_count) VALUES (?,?,?)', (id, expID, sum_counts))
	conn.commit()
	conn.close()

