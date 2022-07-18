#!/usr/bin/env python

import sqlite3
import argparse
import os
import isoSeQL_db as sqlDB
import isoSeQL_parse as fileParse
import sys
import timeit

def main():
	parser=argparse.ArgumentParser(description="parse SQANTI3 output to add isoforms to database")
	parser.add_argument("--classif", help="classification.txt output from SQANTI3")
	parser.add_argument("--genePred", help="genePred from SQANTI3")
	parser.add_argument("--sampleConfig", help="path to sample config file")
	parser.add_argument("--expConfig", help="path to exp config file")
	parser.add_argument("--db", help="path to database")
	parser.add_argument("--sc", default=None, help='path to annotated csv')

	args=parser.parse_args()
	#check for file existence
	errorMsg=""
	if not os.path.isfile(args.classif):
		errorMsg+="\n***EXIT***\n classif file: " + arg.classif + " not found.\n**********\n"
	if not os.path.isfile(args.genePred):
		errorMsg+="\n***EXIT***\n genePred file: " + arg.genePred + " not found.\n**********\n"
	if not os.path.isfile(args.sampleConfig):
		errorMsg+="\n***EXIT***\n sampleConfig file: " + arg.sampleConfig + " not found.\n**********\n"
	if not os.path.isfile(args.expConfig):
		errorMsg+="\n***EXIT***\n expConfig file: " + arg.expConfig + " not found.\n**********\n"
	if args.sc:
		if not os.path.isfile(args.args.sc):
			errorMsg+="\n***EXIT***\n singlecell file: " + arg.sc + " not found.\n**********\n"
	if errorMsg!="":
		print(errorMsg)
		sys.exit(
	else:
		print("\nAll files located.\n")		

	classifInfo=fileParse.parse_classification(args.classif)
	genePredInfo=fileParse.parse_genePred(args.genePred)
	if args.sc:
		scInfo,UMIs=fileParse.parse_singleCell(args.sc)
	else:
		scInfo=None
		UMIs=None
	if not os.path.isfile(args.db):
		sqlDB.make_db(args.db)
	sampleID=sqlDB.addSampleData(args.db, args.sampleConfig)
	expID=sqlDB.addExpData(args.db, args.expConfig, sampleID)
	if expID == "Already in database":
		print("\n***EXIT***\nExperiment and its data have already been added into this database. No isoforms have been added.\n**********\n")
		sys.exit()
	else:
		sqlDB.addIsoforms(args.db, classifInfo, genePredInfo, expID, scInfo, UMIs)





if __name__ =='__main__':
	start=timeit.default_timer()
	main()
	stop=timeit.default_timer()
	print("Complete in {0} sec.".format(stop-start), file=sys.stderr)