#!/usr/bin/env python

import sqlite3
import argparse
import csv
from collections import defaultdict

class Classif:
	def __init__(self,iso_id,cat,gene,transcript,ref_exons,tx_exons,subcat,canonical,count,tx_length):
		self.iso_id = iso_id
		self.cat = cat
		self.gene = gene
		self.transcript = transcript
		self.ref_exons = ref_exons #can be a string if novel so NA
		self.tx_exons = int(tx_exons)
		self.subcat = "variable_ends" if (subcat=="alternative_3end" or subcat=="alternative3end5end" or subcat=="alternative_5end" or subcat=="reference_match") else subcat
		self.canonical = canonical
		self.count = int(count) if (count != "NA") else "NA"
		self.tx_length=int(tx_length)
		self.IEJ = "TRUE" if (self.subcat == "CDS_CDS" or self.subcat == "CDS_UTR" or self.subcat == "UTR_UTR") else "FALSE"
	
	def printObj(self):
		return self.iso_id+'\t'+self.cat+'\t'+self.gene+'\t'+self.transcript+'\t'+str(self.ref_exons)+'\t'+str(self.tx_exons)+'\t'+self.subcat+'\t'+self.canonical+'\t'+str(self.count)+'\t'+self.IEJ+'\t'+str(self.tx_length)+'\n'
		
class Struct:
	def __init__(self,iso_id,chrom,strand,start,end,exStarts,exEnds,exCount):
		self.iso_id = iso_id
		self.chrom = chrom
		self.strand = strand
		self.start = start 
		self.end = end
		self.exStarts = exStarts
		self.exEnds = exEnds
		self.exCount = exCount
		self.exSizes = ','.join([str(self.exEnds[i]-self.exStarts[i]) for i in range(self.exCount)])
		self.junctions = [(self.exEnds[i],self.exStarts[i+1]) for i in range(self.exCount-1)]
		self.exBedStarts=','.join([str(x-self.exStarts[0]) for x in self.exStarts])

	def printObj(self):
		exStart_string = [str(i) for i in self.exStarts]
		exEnd_string = [str(i) for i in self.exEnds]
		return self.iso_id+'\t'+self.chrom+'\t'+self.strand+'\t'+str(self.start)+'\t'+str(self.end)+'\t'+','.join(exStart_string)+'\t'+','.join(exEnd_string)+'\t'+self.exSizes+'\t' + self.exBedStarts+'\n'

# class SingleCell:
# 	def __init__(self, iso_id, barcode, celltype):
# 		self.iso_id = iso_id
# 		self.barcode = barcode
# 		self.celltype = celltype
# 		self.UMIs = set()

# 	def addUMI(self, UMI):
# 		self.UMIs.add(UMI)

# 	def printObj(self):
# 		return self.iso_id+'\t'+self.barcode+'\t'+self.celltype+'\t'+str(self.UMIs)+'\n'

def parse_classification(classif):
	classFile = open(classif, "rt")
	classRead = csv.DictReader(classFile, delimiter="\t")
	classDict={}
	for row in classRead:
		classDict[row['isoform']] = Classif(row['isoform'], row['structural_category'], row['associated_gene'], row['associated_transcript'], row['ref_exons'], row['exons'], row['subcategory'], row['all_canonical'], row['FL'], row['length'])
	return classDict

def parse_genePred(genePred):
	genePredDict={}
	genePredFile = open(genePred, "r")
	while True:
		line=genePredFile.readline().rstrip()
		if not line:
			break
		line=line.split('\t')
		isoID=line[0]
		if isoID in genePredDict.keys():
			print("Repeat transcript " + isoID)
			break
		else: 
			chrom=line[1]
			strand=line[2]
			txStart=int(line[3])
			txEnd=int(line[4])
			exStarts=[int(x) for x in line[8][:-1].split(',')]
			exEnds=[int(x) for x in line[9][:-1].split(',')]
			exCount=int(line[7])
			genePredDict[isoID] = Struct(isoID,chrom,strand,txStart,txEnd,exStarts,exEnds,exCount)
	return genePredDict

def parse_singleCell(sc):
	# scDict={}
	# scFile=open(sc, "r")
	# scRead=csv.DictReader(scFile)
	# for row in scRead:
	# 	isoform=row['pbid']
	# 	barcode=row['BC']
	# 	UMI=row['UMI']
	# 	celltype=row['Celltype']
	# 	if isoform in scDict:
	# 		scDict[isoform][barcode].addUMI(UMI)
	# 	else:
	# 		scDict[isoform]=defaultdict()
	# 		scDict[isoform][barcode]=SingleCell(isoform,barcode,celltype)
	# 		scDict[isoform][barcode].addUMI(UMI)
	# return scDict
	# #want to return two dictionaries - one with just barcode and celltype (to population scInfo table) and then one with sets of UMIs for providing counts
	scFile=open(sc, "r")
	scRead=csv.DictReader(scFile)
	barcodeDict={}
	UMIDict={}
	for row in scRead:
		isoform=row['pbid']
		barcode=row['BC'] #should this be BCrev? b/c that's how cellranger and 10x report it?
		UMI=row['UMI']
		celltype=row['Celltype']
		if isoform in scDict:
			UMIDict[isoform][barcode].add(UMI)
		else:
			UMIDict[isoform]=defaultdict(set)
			UMIDict[isoform][barcode].add(UMI)
		if barcode not in barcodeDict:
			barcodeDict[barcode]=[celltype]
	return barcodeDict,UMIDict

