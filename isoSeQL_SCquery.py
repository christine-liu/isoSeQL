#!/usr/bin/env python

import sqlite3
import csv
import numpy as np
import pandas as pd
import datetime
import argparse
import plotly.graph_objects as go
import plotly.io as pio
import sys
import timeit


def isoprop_plot(db, exp, outPrefix):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	colorDict={"full-splice_match":'#1d2f5f', 'incomplete-splice_match':'#8390FA', 'novel_in_catalog':'#6eaf46', 'novel_not_in_catalog':'#FAC748', 'antisense':'#00bcc0', 'intergenic':'#fa8800', 'genic':'#ac0000', 'genic_intron':'#0096ff', 'fusion':'#cf33ac'}
	#by read
	df_prop=pd.read_sql("SELECT i.category,SUM(c.read_count),s.exp,s.celltype FROM scCounts c INNER JOIN isoform i on i.id=c.isoform_id INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY i.category,s.exp,s.celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals=pd.read_sql("SELECT s.exp, SUM(c.read_count) FROM scCounts c INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY s.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals.rename(columns={'SUM(c.read_count)':'exp_Total'}, inplace=True)
	cell_totals=pd.read_sql("SELECT SUM(c.read_count),s.exp, s.celltype FROM scCounts c INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY s.exp, s.celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	cell_totals.rename(columns={'SUM(c.read_count)':'cell_Total'}, inplace=True)
	prop = df_prop.merge(calc_totals, on=["exp"])
	prop['exp_Proportion']=prop['SUM(c.read_count)']/prop['exp_Total']
	prop=prop.merge(cell_totals,on=["exp","celltype"])
	prop['cell_Proportion']=prop['SUM(c.read_count)']/prop['cell_Total']
	fig = go.Figure()
	pio.full_figure_for_development(fig,warn=False)
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	category2plot=prop.category.unique().tolist()
	colors=[colorDict[i] for i in category2plot]
	for r, col in zip(category2plot, colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.exp_Proportion, name=r, marker_color=col),
		)
	readPlotFile = outPrefix+"_isoCelltypeReadPropPlot.pdf"
	fig.write_image(readPlotFile)
	print("Isoform read proportions per celltype normalized by exp plot saved: " + readPlotFile)
	fig = go.Figure()
	pio.full_figure_for_development(fig,warn=False)
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	for r, col in zip(category2plot, colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.cell_Proportion, name=r, marker_color=col),
		)
	readPlotFile=outPrefix+"_isoCelltypeNormReadPropPlot.pdf"
	fig.write_image(readPlotFile)
	print("Isoform read proportions per celltype normalized by celltype plot saved: " + readPlotFile)
	tableFile=outPrefix+"_isoCelltypeReadsPropTable.txt"
	prop.to_csv(tableFile, sep='\t', index=False, header=True)
	print("Isoform read proportions by celltype table saved: " + tableFile)
#by isoform w/o considering variable ends
	df_prop=pd.read_sql("SELECT category,COUNT(category), exp, celltype FROM (SELECT DISTINCT i.id,i.category, s.exp, s.celltype FROM scCounts c INNER JOIN isoform i on i.id=c.isoform_id INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s))) GROUP BY category,exp,celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	exp_totals=pd.read_sql("SELECT COUNT(isoform_id), exp FROM counts WHERE exp IN (%s) GROUP BY exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	exp_totals.rename(columns={'COUNT(isoform_id)':'exp_Total'}, inplace=True)
	cell_totals=pd.read_sql("SELECT COUNT(isoform_id),exp,celltype FROM (SELECT DISTINCT c.isoform_id, s.exp, s.celltype FROM scCounts c INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s))) GROUP BY exp, celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	cell_totals.rename(columns={'COUNT(isoform_id)':'cell_Total'}, inplace=True)
	prop=df_prop.merge(exp_totals, on=['exp'])
	prop=prop.merge(cell_totals, on=['exp', 'celltype'])
	prop['exp_Proportion']=prop['COUNT(category)']/prop['exp_Total']
	prop['cell_Proportion']=prop['COUNT(category)']/prop['cell_Total']
	fig = go.Figure()
	pio.full_figure_for_development(fig,warn=False)
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	category2plot=prop.category.unique().tolist()
	colors=[colorDict[i] for i in category2plot]
	for r, col in zip(category2plot, colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.exp_Proportion, name=r, marker_color=col),
		)		
	isoPlotFile=outPrefix+"_isoCelltypePropPlot.pdf"
	fig.write_image(isoPlotFile)
	print("Isoform proportions per celltype normalized by exp plot saved: " + readPlotFile)
	fig = go.Figure()
	pio.full_figure_for_development(fig,warn=False)
	fig.update_layout(
		template="simple_white",
		xaxis=dict(title_text="Exp"),
		yaxis=dict(title_text="Proportion"),
		barmode="stack",
	)
	for r, col in zip(category2plot, colors):
		plot_df = prop[prop.category == r]
		fig.add_trace(
			go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.cell_Proportion, name=r, marker_color=col),
		)
	isoPlotFile=outPrefix+"_isoCelltypeNormPropPlot.pdf"
	fig.write_image(isoPlotFile)
	print("Isoform proportions per celltype normalized by celltype plot saved: " + readPlotFile)
	tableFile=outPrefix+"_isoCelltypePropTable.txt"
	prop.to_csv(tableFile, sep='\t', index=False, header=True)
	print("Isoform proportions by celltype table saved: " + tableFile)
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
	df_FSM = pd.read_sql("SELECT tx,gene,SUM(read_count), exp, celltype FROM (SELECT DISTINCT t.tx,t.gene,i.id,c.scID,c.read_count,s.exp,s.celltype FROM scCounts c INNER JOIN isoform i on i.id=c.isoform_id OUTER LEFT JOIN txID t on t.isoform_id=i.id INNER JOIN scInfo s on s.id=c.scID WHERE i.id IN (SELECT id FROM isoform WHERE category=='full-splice_match') AND c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s))) GROUP BY tx,gene,exp,celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	FSM_total=pd.read_sql("SELECT gene,SUM(read_count), exp, celltype FROM (SELECT DISTINCT t.tx,t.gene,i.id,c.scID,c.read_count,s.exp,s.celltype FROM scCounts c INNER JOIN isoform i on i.id=c.isoform_id OUTER LEFT JOIN txID t on t.isoform_id=i.id INNER JOIN scInfo s on s.id=c.scID WHERE i.id IN (SELECT id FROM isoform WHERE category=='full-splice_match') AND c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s))) GROUP BY gene,exp,celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	FSM_total.rename(columns={'SUM(read_count)':'Total'}, inplace=True)
	prop=df_FSM.merge(FSM_total, on=['exp','gene','celltype'])
	prop['Proportion']=prop['SUM(read_count)']/prop['Total']
	for g in gene_list:
		df_gene=prop[(prop["gene"]==g)]
		df_gene.celltype=pd.Categorical(df_gene.celltype, categories=['Ast', 'Ex', 'Inh', 'Mic', 'OPC', 'Oli', 'Per', 'End', 'NA'])
		df_gene=df_gene.sort_values("celltype")
		fig = go.Figure()
		pio.full_figure_for_development(fig,warn=False)
		fig.update_layout(
			template="simple_white",
			xaxis=dict(title_text="Exp"),
			yaxis=dict(title_text="Proportion"),
			barmode="stack",
		)
		for x in df_gene.tx.unique():
			plot_df=df_gene[df_gene.tx==x]
			fig.add_trace(go.Bar(x=[plot_df.exp, plot_df.celltype],y=plot_df.Proportion,name=x))
		fileName=outPrefix+"_FSM_"+g+".pdf"
		fig.write_image(fileName)
		print("FSM read proportions plot saved: " + fileName)
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


def main():
	parser=argparse.ArgumentParser(description="built-in queries to generate plots/tables/visualization of exp/genes of interest")
	subparsers = parser.add_subparsers(dest="subparser_name")
	isoProp_parser = subparsers.add_parser('SCisoProp')
	isoProp_parser.add_argument('--db')
	isoProp_parser.add_argument('--exp')
	isoProp_parser.add_argument('--outPrefix')
	FSM_parser = subparsers.add_parser('SCFSM')
	FSM_parser.add_argument('--db')
	FSM_parser.add_argument('--exp')
	FSM_parser.add_argument('--outPrefix')
	FSM_parser.add_argument('--genes')
	countMat_parser = subparsers.add_parser('SCcountMatrix')
	countMat_parser.add_argument('--db')
	countMat_parser.add_argument('--exp')
	countMat_parser.add_argument('--outPrefix')
	countMat_parser.add_argument('--gene', action='store_true')
	countMat_parser.add_argument('--variable', action='store_true')
	
	args=parser.parse_args()
	if args.subparser_name == "SCisoProp":
		start=timeit.default_timer()
		isoprop_plot(args.db, args.exp, args.outPrefix)
		stop=timeit.default_timer()
		print("Complete in {0} sec.".format(stop-start), file=sys.stderr)
	elif args.subparser_name == "SCFSM":
		start=timeit.default_timer()
		gene_FSM(args.db, args.exp, args.outPrefix, args.genes)
		stop=timeit.default_timer()
		print("Complete in {0} sec.".format(stop-start), file=sys.stderr)
	elif args.subparser_name=="SCcountMatrix":
		start=timeit.default_timer()
		countMatrix(args.db, args.exp, args.outPrefix, args.gene, args.variable)
		stop=timeit.default_timer()
		print("Complete in {0} sec.".format(stop-start), file=sys.stderr)
	else:
		print("\n***ERROR***\n"+args.subparser_name+" is not an option.\n**********\n")




if __name__ =='__main__':
	main()


