#!/usr/bin/env python

import sqlite3
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
import argparse
import plotly.express as px
import plotly.graph_objects as go

def isoprop_plot(db, exp, outPrefix):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	colorDict={"full-splice_match":'#1d2f5f', 'incomplete-splice_match':'#8390FA', 'novel_in_catalog':'#6eaf46', 'novel_not_in_catalog':'#FAC748', 'antisense':'#00bcc0', 'intergenic':'#fa8800', 'genic':'#ac0000', 'genic_intron':'#0096ff', 'fusion':'#cf33ac'}
	#isoform proportion (not counting variable ends)
	df_prop=pd.read_sql("SELECT i.category,COUNT(i.category),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals=pd.read_sql("SELECT COUNT(i.category),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals.rename(columns={'COUNT(i.category)':'Total'}, inplace=True)
	prop=df_prop.merge(calc_totals, on=['exp'])
	prop['Proportion']=prop['COUNT(i.category)']/prop['Total']
	fig = go.Figure()
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	category2plot=prop.category.unique().tolist()
	colors=[colorDict[i] for i in category2plot]
	for r, c in zip(prop.category.unique(), colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=plot_df.exp, y=plot_df.Proportion, name=r, marker_color=c),
    )
	fig.update_layout(xaxis_type='category')
	tableFile=outPrefix+"_isoPropTable.txt"
	prop.to_csv(tableFile, sep='\t', index=True, header=True)
	print("Isoform proportions table saved: " + tableFile)
	plotFile=outPrefix+"_isoPropPlot.pdf"
	fig.write_image(plotFile)
	print("Isoform proportions plot saved: " + plotFile)
	#by read count
	df_exp_read = pd.read_sql("SELECT i.category,SUM(c.read_count),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals = pd.read_sql("SELECT SUM(c.read_count),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals.rename(columns={'SUM(c.read_count)':'Total'}, inplace=True)
	prop=df_exp_read.merge(calc_totals, on=['exp'])
	prop['Proportion']=prop['SUM(c.read_count)']/prop['Total']
	fig = go.Figure()
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	category2plot=prop.category.unique().tolist()
	colors=[colorDict[i] for i in category2plot]
	for r, c in zip(prop.category.unique(), colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=plot_df.exp, y=plot_df.Proportion, name=r, marker_color=c),
    )
	fig.update_layout(xaxis_type='category')
	tableFile=outPrefix+"_isoReadsPropTable.txt"
	prop.to_csv(tableFile, sep='\t', index=True, header=True)
	print("Isoform read proportions table saved: " + tableFile)
	plotFile=outPrefix+"_isoReadsPropPlot.pdf"
	fig.write_image(plotFile)
	print("Isoform read proportions plot saved: " + plotFile)
	return


def make_bed(db, exp, outPrefix, name):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	name_file=open(name, "r")
	name_list=name_file.readlines()
	name_list=[i.rstrip() for i in name_list]
	e=0
	while e < len(exp_list):
		df_bed = pd.read_sql("SELECT DISTINCT i.id,i.chr,x.start,x.end,i.strand,i.category,i.iso_exons,i.junctions, x.ex_sizes,x.ex_starts FROM isoform i INNER JOIN isoform_ends x on i.id=x.isoform_id WHERE x.isoform_id IN (SELECT e.ends_id FROM ends_counts e WHERE e.exp = ?)", conn, params=(exp_list[e],))
		df_bed['name'] = df_bed['id'].map(str)+'_'+df_bed['start'].map(str)+'_'+df_bed['end'].map(str)
		df_bed['score'] = 60
		df_bed['RGB'] = df_bed['category'].apply(lambda x: "29,47,95" if x == "full-splice_match" else '131,144,250' if x=="incomplete-splice_match" else '110,175,70' if x=="novel_in_catalog" else '250,199,72' if x=="novel_not_in_catalog" else '0,188,192' if x=="antisense" else '250,136,0' if x=="intergenic" else '172,0,0' if x=="genic" else '0,150,255' if x=='genic_intron' else '207,51,172' if x=="fusion" else '255,0,0')
		df_bed = df_bed[['chr','start','end','name','score','strand','start','end','RGB','iso_exons','ex_sizes','ex_starts']]
		filename=outPrefix+"_"+name_list[e]+".bed"
		df_bed.to_csv(filename, sep='\t', index=False, header=False)
		header="track visibility=2 itemRgb=On name='"+name_list[e]+"'\n"
		with open(filename, "r+") as f: s=f.read(); f.seek(0); f.write(header + s)
		print("Bedfile saved: " + filename)
		e+=1
	return

def gene_FSM(db, exp, outPrefix, genes):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	gene_file=open(genes, "r")
	gene_list=gene_file.readlines()
	gene_list=[i.rstrip() for i in gene_list]
	df_FSM = pd.read_sql("SELECT t.tx,c.read_count,c.exp,t.gene FROM counts c INNER JOIN txID t on t.isoform_id = c.isoform_id WHERE c.isoform_id IN (SELECT id from isoform WHERE category=='full-splice_match') AND c.exp IN (%s) GROUP BY t.gene,t.tx,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	for g in gene_list:
		df_gene=df_FSM[(df_FSM["gene"]==g)]
		df_gene_pivot=df_gene.pivot(index="exp", columns="tx", values="read_count")
		df_gene_pivot=df_gene_pivot.fillna(0)
		tx_list=df_gene_pivot.columns
		df_gene_pivot['Total'] = df_gene_pivot[tx_list].sum(axis=1)
		pivot_plot=pd.DataFrame()
		for i in tx_list:
			pivot_plot['{}'.format(i)] = df_gene_pivot[i]/df_gene_pivot['Total']	
		ax=pivot_plot.plot.bar(stacked=True).legend(bbox_to_anchor=(1,1), fontsize=8)
		plt.suptitle("Transcript Read Proportion for " + g)
		plt.subplots_adjust(right=0.6)
		fileName=outPrefix+"_FSM_"+g+".pdf"
		plt.savefig(fileName)
		print("FSM read proportions plot saved: " + fileName)
	return

def IEJ_table(db, exp, out):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	df_IEJ=pd.read_sql("SELECT i.id, i.gene, c.exp, c.read_count FROM isoform i INNER JOIN counts c on c.isoform_id = i.id WHERE i.IEJ = 'TRUE' AND c.exp IN (%s) GROUP BY i.gene" % ','.join('?' for i in exp_list), conn, params=exp_list)
	df_IEJ['IEJ_id'] = df_IEJ['gene'] + "_" + df_IEJ['id'].astype(str)
	df_IEJ_pivot=df_IEJ.pivot(index="IEJ_id", columns="exp", values="read_count")
	df_IEJ_pivot=df_IEJ_pivot.fillna(0)
	df_IEJ_pivot.to_csv(out, sep='\t')
	print("IEJ table saved: " + out)
	return

def expInfo(db, out):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	info_table = pd.read_sql("SELECT * FROM exp e INNER JOIN sampleData s on e.sample_id=s.id GROUP BY e.id", conn)
	info_table.columns.values[0]="exp_id"
	info_table = info_table.drop(['id'], axis=1)
	info_table.to_csv(out, sep='\t', index=False, header=True)
	print("Sample info table saved: " + out)
	return

def countMatrix(db, exp, outPrefix, gene=False, variable=False):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	if gene:
		counts = pd.read_sql("SELECT c.isoform_id, c.exp, c.read_count, i.gene FROM counts c INNER JOIN isoform i on i.id=c.isoform_id WHERE exp IN(%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
		geneSum = counts.groupby(['exp', 'gene'])['read_count'].sum().reset_index()
		pivot = geneSum.pivot(index="gene", columns="exp", values="read_count")
		pivot = pivot.fillna(0)
		filename=outPrefix+"_gene_counts_matrix.txt"
		pivot.to_csv(filename, sep='\t', index=True, header=False)
		header="\t"+"\t".join(str(x) for x in pivot.columns.tolist()) + "\n"
		with open (filename, "r+") as f: s=f.read(); f.seek(0); f.write(header+s)
		print("Gene counts matrix saved: " + filename)
		return
	else:
		if not variable:
			commonJxn_counts = pd.read_sql("SELECT isoform_id, exp, read_count FROM counts WHERE exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
			pivot=commonJxn_counts.pivot(index="isoform_id", columns="exp", values="read_count")
			pivot=pivot.fillna(0)
			filename=outPrefix+"_commonJxn_counts_matrix.txt"
			pivot.to_csv(filename, sep='\t', index=True, header=False)
			header="\t"+"\t".join(str(x) for x in pivot.columns.tolist()) + "\n"
			with open(filename, "r+") as f: s=f.read(); f.seek(0); f.write(header + s)
			print("Common Junction counts matrix saved: " + filename)
			return
		else:
			varEnds_counts = pd.read_sql("SELECT x.isoform_id, x.start, x.end, e.exp, e.read_count FROM isoform_ends x INNER JOIN ends_counts e on x.id=e.ends_id WHERE x.id IN (SELECT ends_id FROM ends_counts e WHERE exp IN (%s)) " % ','.join('?' for i in exp_list), conn, params=exp_list)
			varEnds_counts['isoform'] = varEnds_counts['isoform_id'].map(str)+'_'+varEnds_counts['start'].map(str)+'_'+varEnds_counts['end'].map(str)
			pivot = varEnds_counts.pivot(index="isoform", columns="exp", values="read_count")
			pivot=pivot.fillna(0)
			filename=outPrefix+"_variableEnds_counts_matrix.txt"
			pivot.to_csv(filename, sep='\t', index=True, header=False)
			header="\t"+"\t".join(str(x) for x in pivot.columns.tolist()) + "\n"
			with open(filename, "r+") as f: s=f.read(); f.seek(0); f.write(header + s)
			print("Variable ends counts matrix saved: " + filename)
			return

def tappASgff(db, exp, out):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	gff3TX=pd.read_sql("SELECT DISTINCT i.id, i.chr, i.strand, i.gene, i.junctions, i.category, x.start, x.end, x.ex_sizes, t.tx FROM isoform i LEFT OUTER JOIN txID t on i.id=t.isoform_id INNER JOIN isoform_ends x on i.id=x.isoform_id WHERE x.isoform_id IN (SELECT e.ends_id FROM ends_counts e WHERE e.exp IN (%s)) " % ','.join('?' for i in exp_list), conn, params=exp_list)
	gff3TX['isoform'] = gff3TX['id'].map(str)+'_'+gff3TX['start'].map(str)+'_'+gff3TX['end'].map(str)
	gff3TX['tx'].replace(np.nan, "novel", inplace=True)
	outFile = open(out, "w+")
	for i, row in gff3TX.iterrows():
		length=sum(int(x) for x in row['ex_sizes'].split(","))
		txLine=row['isoform']+"\ttappAS\ttranscript\t1\t"+str(length)+"\t.\t"+row['strand']+"\t.\tID="+row['tx']+"; primary_class="+row['category']+'; PosType=T\n'
		a=outFile.write(txLine)
		geneLine=row['isoform']+"\ttappAS\tgene\t1\t"+str(length)+"\t.\t"+row['strand']+"\t.\tID="+row['gene']+"; Name="+row['gene']+"; Desc="+row['gene']+"; PosType=T\n"
		a=outFile.write(geneLine)
		cdsLine=row['isoform']+"\ttappAS\tCDS\t.\t.\t.\t"+row['strand']+"\t.\tID=Protein_"+row['isoform']+"; Name=Protein_"+row['isoform']+'; Desc=Protein_'+row['isoform']+'; PosType=T\n'
		a=outFile.write(cdsLine)
		genomicLine=row['isoform']+"\ttappAS\tgenomic\t1\t1\t.\t"+row['strand']+"\t.\tChr="+row['chr']+'; PosType=G\n'
		a=outFile.write(genomicLine)
		jxn_len=len(row['junctions'])
		if jxn_len > 2:
			jxns=row['junctions'][2:-2].split('), (')
			j=0
			while j < len(jxns):
				startJ=jxns[j].split(", ")[0]
				endJ=jxns[j].split(", ")[1]
				sjLine=row['isoform']+"\ttappAS\tsplice_junction\t"+str(startJ)+"\t"+str(endJ)+"\t.\t"+row['strand']+"\t.\tID_junction_"+str(i)+"; Chr="+row['chr']+"; PosType=G\n"
				a=outFile.write(sjLine)
				j+=1
	print("tappAS-compatible gff3 generated: " + out)
	return	


def main():
	parser=argparse.ArgumentParser(description="built-in queries to generate plots/tables/visualization of exp/genes of interest")
	subparsers = parser.add_subparsers(dest="subparser_name")
	isoProp_parser = subparsers.add_parser('isoProp')
	isoProp_parser.add_argument('--db')
	isoProp_parser.add_argument('--exp')
	isoProp_parser.add_argument('--outPrefix')
	bed_parser = subparsers.add_parser('bed')
	bed_parser.add_argument('--db')
	bed_parser.add_argument('--exp')
	bed_parser.add_argument('--outPrefix')
	bed_parser.add_argument('--name')
	FSM_parser = subparsers.add_parser('FSM')
	FSM_parser.add_argument('--db')
	FSM_parser.add_argument('--exp')
	FSM_parser.add_argument('--outPrefix')
	FSM_parser.add_argument('--genes')
	IEJtab_parser = subparsers.add_parser('IEJtab')
	IEJtab_parser.add_argument('--db')
	IEJtab_parser.add_argument('--exp')
	IEJtab_parser.add_argument('--out')
	sample_parser = subparsers.add_parser('expInfo')
	sample_parser.add_argument('--db')
	sample_parser.add_argument('--out')
	countMat_parser = subparsers.add_parser('countMatrix')
	countMat_parser.add_argument('--db')
	countMat_parser.add_argument('--exp')
	countMat_parser.add_argument('--outPrefix')
	countMat_parser.add_argument('--gene', action='store_true')
	countMat_parser.add_argument('--variable', action='store_true')
	tappASgff_parser = subparsers.add_parser('tappASgff')
	tappASgff_parser.add_argument('--db')
	tappASgff_parser.add_argument('--exp')
	tappASgff_parser.add_argument('--out')

	args=parser.parse_args()
	if args.subparser_name == "isoProp":
		isoprop_plot(args.db, args.exp, args.outPrefix)
	elif args.subparser_name == "bed":
		make_bed(args.db, args.exp, args.outPrefix, args.name)
	elif args.subparser_name == "FSM":
		gene_FSM(args.db, args.exp, args.outPrefix, args.genes)
	elif args.subparser_name == "IEJtab":
		IEJ_table(args.db, args.exp, args.out)
	elif args.subparser_name=="expInfo":
		expInfo(args.db, args.out)
	elif args.subparser_name=="countMatrix":
		countMatrix(args.db, args.exp, args.outPrefix, args.gene, args.variable)
	elif args.subparser_name=="tappASgff":
		tappASgff(args.db, args.exp, args.out)
	else:
		print("\n***ERROR***\n"+args.subparser_name+" is not an option.\n**********\n")




if __name__ =='__main__':
	main()


