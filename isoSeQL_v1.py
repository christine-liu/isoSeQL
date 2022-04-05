#!/usr/bin/env python

#processing of isoseq/cupcake/sqanti data to insert into sqlite database
#inputs: classification.txt, gff, IEJ readnames, a config file with sample information, abundances.txt, database
#parse through files to keep track of isoform information
#insert each isoform into database

#01/28/2021
#used Linnea's 2020jul exosome iso-seq data as test


import sqlite3
import argparse
import csv


class Classif:
	def __init__(self,iso_id,cat,gene,transcript,ref_exons,tx_exons,subcat,canonical,count,tx_length):
		self.iso_id = iso_id
		self.cat = cat
		self.gene = gene
		self.transcript = transcript
		self.ref_exons = ref_exons #can be a string if novel so NA
		self.tx_exons = int(tx_exons)
		self.subcat = subcat
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
		self.exSizes = [self.exEnds[i]-self.exStarts[i] for i in range(self.exCount-1)]
		self.junctions = [(self.exEnds[i],self.exStarts[i+1]) for i in range(self.exCount-1)]

	# def add_exon(self,Exstart,Exend):
	# 	ex_size = int(Exend)- int(Exstart) + 1 #+1 b/c gff doesn't use bed indexing and need to add 1 to get numbrer of bases
	# 	ex_start = int(Exstart) -1 - self.start
	# 	self.Exstarts.append(ex_start)
	# 	self.Exsizes.append(ex_size)
	
	def printObj(self):
		exStart_string = [str(i) for i in self.exStarts]
		exEnd_string = [str(i) for i in self.exEnds]
		exSize_string = [str(i) for i in self.exSizes]
		return self.iso_id+'\t'+self.chrom+'\t'+self.strand+'\t'+str(self.start)+'\t'+str(self.end)+'\t'+','.join(exStart_string)+'\t'+','.join(exEnd_string)+'\t'+','.join(exSize_string)+'\n'

#parse classification text. columns included change depending on options SQANTI3 is run with
def parse_classification(classif):
	print("classif parse")
	# classDict={}
	# classFile = open(classif, "r")
	# colnames = classFile.readline().rstrip().split("\t")
	# iso_index = colnames.index('isoform')
	# cat_index = colnames.index('structural_category')
	# gene_index = colnames.index('associated_gene')
	# tx_index = colnames.index('associated_transcript')
	# refEx_index = colnames.index('ref_exons')
	# txEx_index = colnames.index('exons')
	# subcat_index = colnames.index('subcategory')
	# canonical_index = colnames.index('all_canonical')
	# length_index = colnames.index('length')
	# ###
	# count_index = colnames.index('FL') ##if never chain samples then should just read FL
	# ####
	# while True:
	# 	line = classFile.readline().rstrip()
	# 	if not line:
	# 		break
	# 	line = line.split("\t")
	# 	classDict[line[iso_index]] = Classif(line[iso_index],line[cat_index],line[gene_index],line[tx_index],line[refEx_index],line[txEx_index],line[subcat_index],line[canonical_index],line[count_index],line[length_index])
	# return classDict
	
	#08/11/2021 clean up by using csv dictreader
	classFile = open(classif, "rb")
	classRead = csv.DictReader(classFile, delimiter="\t")
	classDict={}
	for row in classRead:
		classDict[row['isoform']] = Classif(row['isoform'], row['structural_category'], row['associated_gene'], row['associated_transcript'], row['ref_exons'], row['exons'], row['subcategory'], row['all_canonical'], row['FL'], row['length'])
	return classDict

#parse gff output from cDNA_cupcake
def parse_gff(gff):
	print("gff parse")
	gffDict={}
	gffFile = open(gff, "r")
	while True:
		line=gffFile.readline().rstrip()
		if not line:
			break
		line=line.split('\t')
		isoID=line[8].split('transcript_id')[1].split("\"")[1]
		if line[2]=="transcript":
			if isoID in gffDict.keys():
				print "Repeat transcript " + isoID
				break
			else:
				gffDict[isoID] = Struct(isoID,line[0],line[6],line[3],line[4])
		elif line[2]=="exon":
			if isoID not in gffDict.keys():
				print "tx not in dictionary but have exon coords " + isoID
				break
			else:
				gffDict[isoID].add_exon(line[3],line[4])
		else:
			continue
	return gffDict

#parse genePred file from SQANTI
def parse_genePred(genePred):
	print("genePred parse")
	genePredDict={}
	genePredFile = open(genePred, "r")
	while True:
		line=genePredFile.readline().rstrip()
		if not line:
			break
		line=line.split('\t')
		isoID=line[0]
		if isoID in genePredDict.keys():
			print "Repeat transcript " + isoID
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

# #update IEJ status
# def parse_IEJ(IEJ, classifObj):
# 	IEJFile = open(IEJ, "r")
# 	while True:
# 		line=IEJFile.readline().rstrip()
# 		if not line:
# 			break
# 		classifObj[line].IEJ=True
# 	return
##07/14/2021 no longer have separate IEJ analysis that outputs like this


def main():
	parser=argparse.ArgumentParser(description="parse SQANTI3 output to add isoforms to database")
	parser.add_argument("--classif", help="classification.txt output from SQANTI3")
	# parser.add_argument("--gff", help="gff output from SQANTI3")
	parser.add_argument("--genePred", help="genePred from SQANTI3")
	# parser.add_argument("--IEJ", help="IEJ readnames (isoforms) from IEJ script")
	#parser.add_argument("--config", help="config file with sample information")
	#parser.add_argument("--db", help="path to database to populate with isoform run info")
	# print("HI") 
	args=parser.parse_args()
	classifInfo=parse_classification(args.classif) #dictionary with classification objects
	#gffInfo=parse_gff(args.gff) #dictionary with structure info objects
	genePredInfo=parse_genePred(args.genePred)
	# parse_IEJ(args.IEJ,classifInfo)

	# #for testing purposes, print dictionaries to file to see if information is correct?
	# classifPrint=open("testClassifDict.txt", "w+")
	classifPrint=open("2021aug11_classifTestDict.txt", "w+")
	for key in classifInfo.keys():
		classifPrint.write(classifInfo[key].printObj())

	# gffPrint=open("testGFFDict.txt", "w+")
	# for key in gffInfo.keys():
	# 	gffPrint.write(gffInfo[key].printObj())
	genePredPrint=open("2021aug11_genePredTestDict.txt", "w+")
	for key in genePredInfo.keys():
		genePredPrint.write(genePredInfo[key].printObj())

	# print("keys should be the same...")
	# print(classifInfo.keys().sort()==gffInfo.keys().sort())
	# print("finished")

if __name__ =='__main__':
	main()



