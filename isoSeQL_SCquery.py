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
import seaborn as sns
import plotly.graph_objects as go
import plotly.io as pio
pio.full_figure_for_development(fig,warn=False)

def isoprop_plot(db, exp, outPrefix):
	conn=sqlite3.connect(db)
	c=conn.cursor()
	exp_file=open(exp, "r")
	exp_list=exp_file.readlines()
	exp_list=[i.rstrip() for i in exp_list]
	#by read
	df_prop=pd.read_sql("SELECT i.category,SUM(c.read_count),s.exp,s.celltype FROM scCounts c LEFT OUTER JOIN isoform i on i.id=c.isoform_id INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY i.category,s.exp,s.celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals=pd.read_sql("SELECT s.exp, SUM(c.read_count) FROM scCounts c INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY s.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	calc_totals.rename(columns={'SUM(c.read_count)':'Total'}, inplace=True)
	prop = df_prop.merge(calc_totals, on=["exp"])
	prop['Proportion']=prop['SUM(c.read_count)']/prop['Total']
	prop['exp_celltype']=prop['exp'].map(str)+'_'+prop['celltype']
	cell_totals=pd.read_sql("SELECT SUM(c.read_count),s.exp, s.celltype FROM scCounts c INNER JOIN scInfo s on s.id=c.scID WHERE c.scID IN (SELECT id FROM scInfo WHERE exp IN (%s)) GROUP BY s.exp, s.celltype" % ','.join('?' for i in exp_list), conn, params=exp_list)
	cell_totals['exp_celltype']=cell_totals['exp'].map(str)+'_'+cell_totals['celltype']
	cell_totals.rename(columns={'SUM(c.read_count)':'cellTotal'}, inplace=True)
	cell_totals=cell_totals[["cellTotal","exp_celltype"]]
	prop=prop.merge(cell_totals,on="exp_celltype")
	prop['celltypeProportion']=prop['SUM(c.read_count)']/prop['cellTotal']
	fig = go.Figure()
	fig.update_layout(
	    template="simple_white",
	    xaxis=dict(title_text="Exp"),
	    yaxis=dict(title_text="Proportion"),
	    barmode="stack",
	)
	colorDict={"full-splice_match":'#1d2f5f', 'incomplete-splice_match':'#8390FA', 'novel_in_catalog':'#6eaf46', 'novel_not_in_catalog':'#FAC748', 'antisense':'#00bcc0', 'intergenic':'#fa8800', 'genic':'ac0000', 'genic_intron':'#0096ff', 'fusion':'cf33ac'}
	category2plot=prop.category.unique().tolist()
	colors=[colorDict[i] for i in category2plot]
	for r, c in zip(category2plot, colors):
	    plot_df = prop[prop.category == r]
	    fig.add_trace(
	        go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.Proportion, name=r, marker_color=c),
	    )
	readPlotFile = outPrefix+"_isoCelltypeReadPropPlot.pdf"
	fig.write_image(readPlotFile)
	print("Isoform read proportions plot saved: " + readPlotFile)
	fig = go.Figure()
	fig.update_layout(
	    template="simple_white",
	    xaxis=dict(title_text="Exp"),
	    yaxis=dict(title_text="Proportion"),
	    barmode="stack",
	)
	for r, c in zip(category2plot, colors):
	    plot_df = prop[prop.category == r]
	    fig.add_trace(
	        go.Bar(x=[plot_df.exp, plot_df.celltype], y=plot_df.celltypeProportion, name=r, marker_color=c),
	    )
	readPlotFile=outPrefix+"_isoCelltypeNormReadPropPlot.pdf"
	fig.write_image(readPlotFile)
	print("Isoform read proportions plot saved: " + readPlotFile)
	df_bed['name'] = df_bed['id'].map(str)+'_'+df_bed['start'].map(str)+'_'+df_bed['end'].map(str)




################

	df_prop=pd.read_sql("SELECT i.category,COUNT(i.category),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	pivot=df_prop.pivot(index="exp", columns="category", values="COUNT(i.category)")
	cats = pivot.columns.tolist()
	# cats = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'antisense', 'genic', 'intergenic', 'fusion', 'genic_intron']
	pivot['Total'] = pivot[cats].sum(axis=1)
	pivot_plot=pd.DataFrame()
	for i in cats:
		pivot_plot['{}'.format(i)] = pivot[i]/pivot['Total']
	tableFile=outPrefix+"_isoPropTable.txt"
	pivot_plot.to_csv(tableFile, sep='\t', index=True, header=True)
	print("Isoform proportions table saved: " + tableFile)
	plotFile=outPrefix+"_isoPropPlot.pdf"
	ax=pivot_plot.plot.bar(stacked=True, color={'full-splice_match':'#006e00', 'incomplete-splice_match':'#008cf9', 'novel_in_catalog':'#b80058', 'novel_not_in_catalog':'#ebac23', 'antisense':'#00bbad', 'genic':'#d163e6', 'intergenic':'#b24502', 'fusion':'#ff9287', 'genic_intron':'#5954d6'}).legend(bbox_to_anchor=(1,1), fontsize=8)
	plt.subplots_adjust(right=0.6)
	plt.savefig(plotFile)
	print("Isoform proportions plot saved: " + plotFile)
	df_exp_read = pd.read_sql("SELECT i.category,SUM(c.read_count),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
	pivot_read = df_exp_read.pivot(index="exp", columns="category", values="SUM(c.read_count)")
	cats_reads = pivot_read.columns.tolist()
	pivot_read['Total'] = pivot_read[cats_reads].sum(axis=1)
	pivot_filterReads=pd.DataFrame()
	for i in cats:
		pivot_filterReads['{}'.format(i)] = pivot_read[i]/pivot_read['Total']
	readTableFile=outPrefix+"_isoReadsPropTable.txt"
	pivot_filterReads.to_csv(readTableFile, sep='\t', index=True, header=True)
	print("Isoform read proportions table saved: " + readTableFile)
	readPlotFile = outPrefix+"_isoReadPropPlot.pdf"
	ax=pivot_filterReads.plot.bar(stacked=True, color={'full-splice_match':'#006e00', 'incomplete-splice_match':'#008cf9', 'novel_in_catalog':'#b80058', 'novel_not_in_catalog':'#ebac23', 'antisense':'#00bbad', 'genic':'#d163e6', 'intergenic':'#b24502', 'fusion':'#ff9287', 'genic_intron':'#5954d6'}).legend(bbox_to_anchor=(1,1), fontsize=8)
	plt.subplots_adjust(right=0.6)
	plt.savefig(readPlotFile)
	print("Isoform read proportions plot saved: " + readPlotFile)
	return


###
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
	isoProp_parser = subparsers.add_parser('isoProp')
	isoProp_parser.add_argument('--db')
	isoProp_parser.add_argument('--exp')
	isoProp_parser.add_argument('--outPrefix')
	FSM_parser = subparsers.add_parser('FSM')
	FSM_parser.add_argument('--db')
	FSM_parser.add_argument('--exp')
	FSM_parser.add_argument('--outPrefix')
	FSM_parser.add_argument('--genes')
	countMat_parser = subparsers.add_parser('countMatrix')
	countMat_parser.add_argument('--db')
	countMat_parser.add_argument('--exp')
	countMat_parser.add_argument('--outPrefix')
	countMat_parser.add_argument('--gene', action='store_true')
	countMat_parser.add_argument('--variable', action='store_true')
	
	args=parser.parse_args()
	if args.subparser_name == "isoProp":
		isoprop_plot(args.db, args.exp, args.outPrefix)
	elif args.subparser_name == "FSM":
		gene_FSM(args.db, args.exp, args.outPrefix, args.genes)
	elif args.subparser_name=="countMatrix":
		countMatrix(args.db, args.exp, args.outPrefix, args.gene, args.variable)
	else:
		print("\n***ERROR***\n"+args.subparser_name+" is not an option.\n**********\n")




if __name__ =='__main__':
	main()


