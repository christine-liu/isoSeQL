#!/usr/bin/env python

import sqlite3
import csv

def make_db(database):
	conn=sqlite3.connect(database)
	c=conn.cursor()
	c.execute("CREATE TABLE versionInfo (name TEXT PRIMARY KEY, visoSeQL TEXT)")
	c.execute("INSERT INTO versionInfo(name, visoSeQL) VALUES (?,?)", (database, 'v1.0.1'))
	c.execute("CREATE TABLE isoform (id INTEGER PRIMARY KEY, chr TEXT, strand TEXT, junctions TEXT, gene TEXT, iso_exons INTEGER, subcategory TEXT, canonical TEXT, IEJ INTEGER, category TEXT, UNIQUE(chr, strand, junctions, gene))")
	c.execute("CREATE TABLE isoform_ends (id INTEGER PRIMARY KEY, isoform_id INTEGER, chr TEXT, start INTEGER, end INTEGER, ex_sizes TEXT, ex_starts TEXT, FOREIGN KEY(isoform_id) REFERENCES isoform(id), UNIQUE(isoform_id, chr, start, end))")
	c.execute("CREATE TABLE counts (id INTEGER PRIMARY KEY, isoform_id INTEGER, exp INTEGER, read_count INTEGER, FOREIGN KEY(isoform_id) REFERENCES isoform(id),FOREIGN KEY(exp) REFERENCES exp(id), UNIQUE(isoform_id, exp))")
	c.execute("CREATE TABLE ends_counts (ends_id INTEGER, exp INTEGER, read_count INTEGER, FOREIGN KEY(ends_id) REFERENCES isoform_ends(id), FOREIGN KEY(exp) REFERENCES exp(id), UNIQUE(ends_id, exp, read_count))")
	c.execute("CREATE TABLE exp (id INTEGER PRIMARY KEY, sample_id INTEGER, RIN REAL, seq_date DATE, platform TEXT, method TEXT, vMap TEXT, vReference TEXT, vAnnot TEXT, vLima TEXT, vCCS TEXT, vIsoseq3 TEXT, vCupcake TEXT, vSQANTI TEXT, exp_name TEXT, FOREIGN KEY(sample_id) REFERENCES sampleData(id), UNIQUE(sample_id, RIN, seq_date, platform, method, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI, exp_name))")
	c.execute("CREATE TABLE sampleData (id INTEGER PRIMARY KEY, sample_name TEXT, tissue TEXT, disease TEXT, age INTEGER, sex TEXT,UNIQUE(sample_name, tissue, disease, age, sex))")
	c.execute("CREATE TABLE PBID (id INTEGER PRIMARY KEY, PBID TEXT, exp INTEGER, isoform_id INTEGER, FOREIGN KEY(exp) REFERENCES exp(id), FOREIGN KEY(isoform_id) REFERENCES isoform(id), UNIQUE(PBID, exp, isoform_id))")
	c.execute("CREATE TABLE txID (id INTEGER PRIMARY KEY, tx TEXT, exp INTEGER, isoform_id INTEGER, gene TEXT, FOREIGN KEY(exp) REFERENCES exp(id), FOREIGN KEY(isoform_id) REFERENCES isoform(id), UNIQUE(tx, exp, isoform_id))")
	c.execute("CREATE TABLE scInfo (id INTEGER PRIMARY KEY, exp INTEGER, barcode TEXT, celltype TEXT, FOREIGN KEY(exp) REFERENCES exp(id), UNIQUE(id, exp, barcode, celltype))")
	c.execute("CREATE TABLE scCounts (isoform_id INTEGER, scID INTEGER, read_count INTEGER, FOREIGN KEY(isoform_id) REFERENCES isoform(id), FOREIGN KEY(scID) REFERENCES scInfo(id), UNIQUE(isoform_id, scID, read_count))")
	c.execute("CREATE TABLE scCounts_ends (ends_id INTEGER, scID INTEGER, read_count INTEGER, FOREIGN KEY(ends_id) REFERENCES isoform_ends(id), FOREIGN KEY(scID) REFERENCES scInfo(id), UNIQUE(ends_id, scID, read_count))")

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
		c.execute('SELECT id FROM exp WHERE sample_id = ? AND RIN = ? AND seq_date = ? AND platform = ? AND method = ? AND vMap = ? AND vReference = ? AND vAnnot = ? AND vLima = ? AND vCCS = ? AND vIsoseq3 = ? AND vCupcake = ? AND vSQANTI = ? AND exp_name = ?', (sampleID, row['RIN'], row['date'], row['platform'], row['method'], row['vMap'], row['vReference'], row['vAnnot'], row['vLima'], row['vCCS'], row['vIsoseq3'], row['vCupcake'], row['vSQANTI'], row['exp_name'],))
		expID = c.fetchall()
		if len(expID) == 0:
			c.execute("INSERT INTO exp(sample_id, RIN, seq_date, platform, method, vMap, vReference, vAnnot, vLima, vCCS, vIsoseq3, vCupcake, vSQANTI, exp_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (sampleID, row['RIN'], row['date'], row['platform'], row['method'], row['vMap'], row['vReference'], row['vAnnot'], row['vLima'], row['vCCS'], row['vIsoseq3'], row['vCupcake'], row['vSQANTI'], row['exp_name'],))
			expID=c.lastrowid
		else:
			expID = "Already in database"
	conn.commit()
	conn.close()
	return expID

def addIsoforms(database, classif, genePred, expID, scInfo=None, UMIs=None):
	conn=sqlite3.connect(database)
	c=conn.cursor()
	observedIso=set()
	if scInfo: #add all barcode/celltype info into scInfo table
			for barcode in scInfo.keys():
				c.execute('SELECT id FROM scInfo WHERE exp = ? AND barcode = ? AND celltype = ?', (expID, barcode, scInfo[barcode]))
				scID = c.fetchall()
				if len(scID) == 0:
					c.execute('INSERT INTO scInfo(exp, barcode, celltype) VALUES (?,?,?)', (expID, barcode, scInfo[barcode],))
	for iso in classif.keys():
		c.execute('SELECT id FROM isoform WHERE chr = ? AND strand = ? AND junctions = ? and gene = ?', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions), classif[iso].gene)) #check if isoform is already in table 
		isoID = c.fetchall()
		if len(isoID) == 0:
			c.execute('INSERT INTO isoform(chr, strand, junctions, gene, iso_exons, subcategory, canonical, IEJ, category) VALUES (?,?,?,?,?,?,?,?,?)', (genePred[iso].chrom, genePred[iso].strand, str(genePred[iso].junctions), classif[iso].gene, classif[iso].tx_exons, classif[iso].subcat, classif[iso].canonical, classif[iso].IEJ, classif[iso].cat,)) #if not add into table and return id 
			isoID=c.lastrowid
		else:
			isoID=isoID[0][0] #if so, return id
		observedIso.add(isoID)
		c.execute('SELECT id FROM isoform_ends WHERE isoform_id = ? AND chr = ? AND start = ? AND end = ? AND ex_sizes = ? AND ex_starts = ?', (isoID, genePred[iso].chrom, genePred[iso].start, genePred[iso].end, genePred[iso].exSizes, genePred[iso].exBedStarts)) #check if exact isoform already in table
		isoEndID = c.fetchall()
		if len(isoEndID) == 0:
			c.execute('INSERT INTO isoform_ends(isoform_id, chr, start, end, ex_sizes, ex_starts) VALUES (?,?,?,?,?,?)', (isoID, genePred[iso].chrom, genePred[iso].start, genePred[iso].end, genePred[iso].exSizes, genePred[iso].exBedStarts)) #if not add into table and return id
			isoEndID=c.lastrowid
		else:
			isoEndID=isoEndID[0][0] #if so, return id
		c.execute('INSERT INTO ends_counts(ends_id, exp, read_count) VALUES (?,?,?)', (isoEndID, expID, classif[iso].count)) #should be no repeats of exact isoforms in an exp, so can just directly add counts to table
		c.execute('INSERT INTO PBID(PBID, exp, isoform_id) VALUES (?,?,?)', (iso, expID, isoID)) #add in PBID info to match back to original SQANTI3 output files
		if classif[iso].cat == "full-splice_match" or classif[iso].cat=="incomplete-splice_match":
			c.execute('INSERT OR IGNORE INTO txID(isoform_id, exp, tx, gene) VALUES (?,?,?,?)', (isoID, expID, classif[iso].transcript, classif[iso].gene)) #keep track of txIDs
		if UMIs:
			if iso in UMIs.keys():
				for barcode in UMIs[iso].keys():
					c.execute('SELECT id FROM scInfo WHERE exp = ? AND barcode = ?', (expID, barcode))
					scID = c.fetchall()
					if len(scID) == 0:
						print("ERROR: missing single cell info")
						exit
					else:
						scID = scID[0][0]
						c.execute('INSERT INTO scCounts_ends(ends_id, scID, read_count) VALUES (?,?,?)', (isoEndID, scID, len(UMIs[iso][barcode],)))
						#need to check then for each isoform and where to add the counts, all isoforms should be in the database b/c parse classif first

						#is it any faster to parse through all isoforms and then each cell associated with each isoform vs each cell and then each isoform for every cell?
			else:
				continue


	for isoform in observedIso:
		#query counts from matching isoform_ends ids and sum to get counts to add to counts table (and correspondingly for single-cell info)
		c.execute('SELECT SUM(read_count) FROM ends_counts WHERE exp = ? AND ends_id IN (SELECT id FROM isoform_ends WHERE isoform_id = ?)', (expID, isoform,))
		sum_counts = int(c.fetchall()[0][0])
		c.execute('INSERT INTO counts(isoform_id, exp, read_count) VALUES (?,?,?)', (isoform, expID, sum_counts,)) #now that have counted all exact versions of isoform, can add counts
		if UMIs:
			c.execute('SELECT scID,SUM(read_count) FROM scCounts_ends WHERE scID IN (SELECT id FROM scInfo WHERE exp = ?) AND ends_id IN (SELECT id FROM isoform_ends WHERE isoform_id = ?) GROUP BY scID', (expID, isoform,))
			count_list=c.fetchall()
			for cell in count_list:
				c.execute("INSERT INTO scCounts(isoform_id, scID, read_count) VALUES (?,?,?)", (isoform, cell[0], cell[1],))




	conn.commit()
	conn.close()

