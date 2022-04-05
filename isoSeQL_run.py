#!/usr/bin/env python

import sqlite3
import argparse
import os
import isoSeQL_db as sqlDB
import isoSeQL_parse as fileParse
import sys


def main():
	parser=argparse.ArgumentParser(description="parse SQANTI3 output to add isoforms to database")
	parser.add_argument("--classif", help="classification.txt output from SQANTI3")
	parser.add_argument("--genePred", help="genePred from SQANTI3")
	parser.add_argument("--sampleConfig", help="path to sample config file")
	parser.add_argument("--expConfig", help="path to exp config file")
	parser.add_argument("--db", help="path to database")

	args=parser.parse_args()
	classifInfo=fileParse.parse_classification(args.classif)
	genePredInfo=fileParse.parse_genePred(args.genePred)

	if not os.path.isfile(args.db):
		sqlDB.make_db(args.db)
	sampleID=sqlDB.addSampleData(args.db, args.sampleConfig)
	expID=sqlDB.addExpData(args.db, args.expConfig, sampleID)
	if expID == "Already in database":
		print("\n***EXIT***\nExperiment and its data have already been added into this database. No isoforms have been added.\n**********\n")
		sys.exit()
	else:
		sqlDB.addIsoforms(args.db, classifInfo, genePredInfo, expID)





if __name__ =='__main__':
	main()