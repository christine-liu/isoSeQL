import sqlite3
import csv

conn=sqlite3.connect("20210901_test.db")
c=conn.cursor()

columns=[column[0] for column in c.description]
results=[]
c.execute("SELECT * from exp")
for row in c.fetchall():
	results.append(dict(zip(columns, row)))

with open("testExp.csv", "w") as new_file:
	fieldnames=columns
	writer=csv.DictWriter(new_file, fieldnames=fieldnames, delimiter='\t')
	writer.writeheader()
	for line in results:
		writer.writerow(line)

columns=[column[0] for column in c.description]
results=[]
c.execute("SELECT * from sampleData")
for row in c.fetchall():
	results.append(dict(zip(columns, row)))

with open("testSampleData.csv", "w") as new_file:
	fieldnames=columns
	writer=csv.DictWriter(new_file, fieldnames=fieldnames, delimiter='\t')
	writer.writeheader()
	for line in results:
		writer.writerow(line)

c.execute("SELECT e.id,s.patient,s.tissue,s.disease,s.RIN,s.age,e.seq_date FROM sampleData s, exp e WHERE s.id = e.patient_id GROUP BY s.id")
columns=[column[0] for column in c.description]
results=[]
for row in c.fetchall():
	results.append(dict(zip(columns, row)))

with open("testSampleDataExp.csv", "w") as new_file:
	fieldnames=columns
	writer=csv.DictWriter(new_file, fieldnames=fieldnames, delimiter='\t')
	writer.writeheader()
	for line in results:
		writer.writerow(line)

#09/20/2021 making a plot of proportions of isoforms in different categories (by read or by isoform? both?)
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

conn = sqlite3.connect("20210901_test.db")
c=conn.cursor()
df = pd.read_sql_query("SELECT category,COUNT(category) FROM isoform WHERE category IN ('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'antisense', 'genic', 'intergenic', 'fusion', 'genic_intron') GROUP BY category", conn)
cat_list=['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'antisense', 'genic', 'intergenic', 'fusion', 'genic_intron']
df = pd.read_sql("SELECT category,COUNT(category) FROM isoform WHERE category IN (%s) GROUP BY category" % ','.join('?' for i in cat_list), conn, params=cat_list)

df = pd.read_sql("SELECT category,COUNT(category) FROM isoform WHERE category IN (%s) GROUP BY category" % ','.join('?' for i in cat_list), conn, params=cat_list)
exp_list=[1,2]
# df_exp = pd.read_sql("SELECT i.category,COUNT(i.category),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) AND i.category IN (?,?,?,?,?,?,?,?,?) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list+cat_list)
df_exp = pd.read_sql("SELECT i.category,COUNT(i.category),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
pivot=df_exp.pivot(index="exp", columns="category", values="COUNT(i.category)")
cats = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'antisense', 'genic', 'intergenic', 'fusion', 'genic_intron']
pivot['Total'] = pivot[cats].sum(axis=1)
for i in cats:
	pivot['{}_Percent'.format(i)] = pivot[i]/pivot['Total']
pivot_filter=pivot[["full-splice_match_Percent", "incomplete-splice_match_Percent", "novel_in_catalog_Percent", "novel_not_in_catalog_Percent"]]
# ax=pivot_filter.plot.bar(stacked=True, color=['#1D2F6F', '#8390FA', '#6EAF46', '#FAC748']).legend(bbox_to_anchor=(1,1), fontsize=8)
ax=pivot_filter.plot.bar(stacked=True, color={'full-splice_match_Percent':'#006e00', 'incomplete-splice_match_Percent':'#008cf9', 'novel_in_catalog_Percent':'#b80058', 'novel_not_in_catalog_Percent':'#ebac23', 'antisense_Percent':'#00bbad', 'genic_Percent':'#d163e6', 'intergenic_Percent':'#b24502', 'fusion_Percent':'#ff9287', 'genic_intron_Percent':'#5954d6'}).legend(bbox_to_anchor=(1,1), fontsize=8)
plt.subplots_adjust(right=0.60)
plt.savefig("/mnt/data/chunlab/Christine/isoSeQL_TEST/20210921_isoformcat_noEnds.pdf")

df_exp_read = pd.read_sql("SELECT i.category,SUM(c.read_count),c.exp FROM counts c INNER JOIN isoform i on i.id = c.isoform_id WHERE c.exp IN (%s) GROUP BY i.category,c.exp" % ','.join('?' for i in exp_list), conn, params=exp_list)
pivot_read = df_exp_read.pivot(index="exp", columns="category", values="SUM(c.read_count)")
cats = ['full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog']
pivot_read['Total'] = pivot_read[cats].sum(axis=1)
pivot_filterReads=pd.DataFrame()
for i in cats:
	pivot_filterReads['{}'.format(i)] = pivot_read[i]/pivot_read['Total']
ax=pivot_filterReads.plot.bar(stacked=True, color={'full-splice_match':'#006e00', 'incomplete-splice_match':'#008cf9', 'novel_in_catalog':'#b80058', 'novel_not_in_catalog':'#ebac23', 'antisense':'#00bbad', 'genic':'#d163e6', 'intergenic':'#b24502', 'fusion':'#ff9287', 'genic_intron':'#5954d6'}).legend(bbox_to_anchor=(1,1), fontsize=8)
plt.subplots_adjust(right=0.6)
plt.savefig("/mnt/data/chunlab/Christine/isoSeQL_TEST/20210921_isoformcat_reads.pdf")


ax=pivot_filterReads.plot.bar(stacked=True, color={'full-splice_match':'#006e00', 'incomplete-splice_match':'#008cf9', 'novel_in_catalog':'#B24502', 'novel_not_in_catalog':'#ebac23', 'antisense':'#00bbad', 'genic':'#d163e6', 'intergenic':'#b24502', 'fusion':'#ff9287', 'genic_intron':'#5954d6'}).legend(bbox_to_anchor=(1,1), fontsize=8)


#https://towardsdatascience.com/stacked-bar-charts-with-pythons-matplotlib-f4020e4eb4a7
# def plot_stackedbar_p(df, labels, colors, title, subtitle):
# 	fig=plt.figure()
# 	fields = df.columns.tolist()
# 	#figure and axis
# 	fig, ax = plt.subplots(1, figsize=(12,10))
# 	#plot bars
# 	left = len(df) * [0]
# 	for idx, name in enumerate(fields):
# 		plt.barh(df.index, df[name], left=left, color=colors[idx])
# 		left=left+df[name]
# 	#title and subtitle
# 	plt.title(title, loc='left')
# 	plt.text(0, ax.get_yticks()[-1]+ 0.75, subtitle)
# 	#legend
# 	plt.legend(labels, bbox_to_anchor=([0.58,1,0,0]),ncol=4,frameon=False)
# 	#remove spines
# 	ax.spines['right'].set_visible(False)
# 	ax.spines['left'].set_visible(False)
# 	ax.spines['top'].set_visible(False)
# 	ax.spines['bottom'].set_visible(False)
# 	#format x ticks
# 	xticks=np.arange(0,1.1,0.1)
# 	xlabels=['{}%'.format(i) for i in np.arange(0,101,10)]
# 	plt.xticks(xticks,xlabels)
# 	#adjust limits and draw grid lines
# 	plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
# 	ax.xaxis.grid(color="gray", linestyle="dashed")
# 	plt.close(fig)
# 	plt.savefig('/mnt/data/chunlab/Christine/isoSeQL_TEST/20210920_isoformcat_noEnds.pdf')

# labels=['full-splice_match', 'incomplete-splice_match','novel_in_catalog', 'novel_not_in_catalog']
# colors=['#1D2F6F', '#8390FA', '#6EAF46', '#FAC748']
# title="Isoforms by Sample and Structural Category"
# subtitle="Proportion of Total Isoforms by Structural Category"
# plot_stackedbar_p(pivot_filter, labels, colors, title, subtitle)

#09/21/2021 
#query info to make a bed file
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

conn = sqlite3.connect("20210923_test.db")
c=conn.cursor()

exp_list=[1]
df_bed = pd.read_sql("SELECT i.id,i.chr,x.start,x.end,i.strand,i.category,i.iso_exons,i.junctions, x.ex_sizes, x.ex_starts FROM isoform i INNER JOIN isoform_ends x on i.id = x.isoform_id WHERE x.exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
df_bed['name'] = df_bed['id'].map(str)+'_'+df_bed['start'].map(str)+'_'+df_bed['end'].map(str)
df_bed['score'] = 60
df_bed['RGB'] = df_bed['category'].apply(lambda x: "29,47,95" if x == "full-splice_match" else '131,144,250' if x=="incomplete-splice_match" else '110,175,70' if x=="novel_in_catalog" else '250,199,72' if x=="novel_not_in_catalog" else '0,188,192' if x=="antisense" else '250,136,0' if x=="intergenic" else '172,0,0' if x=="genic" else '0,150,255' if x=='genic_intron' else '207,51,172' if x=="fusion" else '255,0,0')
df_bed = df_bed[['chr','start','end','name','score','strand','start','end','RGB','iso_exons','ex_sizes','ex_starts']]
df_bed.to_csv("/mnt/data/chunlab/Christine/isoSeQL_TEST/20210923_bedTest.bed", sep='\t', index=False, header=False)



#try with one of the ALS datasets
#2714
conn = sqlite3.connect("ALS_2714test/20210923_ALS2714_test.db")
c=conn.cursor()

exp_list=[1]
df_bed = pd.read_sql("SELECT DISTINCT i.id,i.chr,x.start,x.end,i.strand,i.category,i.iso_exons,i.junctions, x.ex_sizes, x.ex_starts FROM isoform i INNER JOIN isoform_ends x on i.id = x.isoform_id WHERE x.exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
#if give multiple exp and an isoform is in multiple exp then only draw once. (DISTINCT keyword). doesn't check which ones are in common or unique to particular exp
df_bed['name'] = df_bed['id'].map(str)+'_'+df_bed['start'].map(str)+'_'+df_bed['end'].map(str)
df_bed['score'] = 60
df_bed['RGB'] = df_bed['category'].apply(lambda x: "29,47,95" if x == "full-splice_match" else '131,144,250' if x=="incomplete-splice_match" else '110,175,70' if x=="novel_in_catalog" else '250,199,72' if x=="novel_not_in_catalog" else '0,188,192' if x=="antisense" else '250,136,0' if x=="intergenic" else '172,0,0' if x=="genic" else '0,150,255' if x=='genic_intron' else '207,51,172' if x=="fusion" else '255,0,0')
df_bed = df_bed[['chr','start','end','name','score','strand','start','end','RGB','iso_exons','ex_sizes','ex_starts']]
df_bed.to_csv("/mnt/data/chunlab/Christine/isoSeQL_TEST/ALS_2714test/20210923_ALS2714_test.bed", sep='\t', index=False, header=False)
with open("20210923_ALS2714_test.bed", "r+") as f: s=f.read(); f.seek(0); f.write("track visibility=2 itemRgb=On name='ALS2714'\n" + s)

#idk try select from list of exp and list of genes? but if multiple genes, then separate plots, so maybe that's just a loop??
#in ALS_2714_C030
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

conn=sqlite3.connect("20211008_2ALSsamp.db")
c=conn.cursor()

gene_list = ["MAPK10", "BIN1", "NDRG4"]
exp_list=[1,2]
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
	fileName="/mnt/data/chunlab/Christine/isoSeQL_TEST/ALS_2714_C030/20211008_FSM_"+g+".pdf"
	plt.savefig(fileName)

#generate table of isoforms with IEJs, which samples they appeared in and how many reads
#does not check if the ends of the isoforms are exactly the same but checks if all the jxns are
#isoforms that contain the same IEJ but vary with other junctions will not be counted as the same ...I think this still makes sense b/c we're not looking at the specific junction but the isoform overall
conn=sqlite3.connect("20211008_2ALSsamp.db")
c=conn.cursor()

exp_list=[1,2]
df_IEJ=pd.read_sql("SELECT i.id, i.gene, c.exp, c.read_count FROM isoform i INNER JOIN counts c on c.isoform_id = i.id WHERE i.IEJ = 'TRUE' AND c.exp IN (%s) GROUP BY i.gene" % ','.join('?' for i in exp_list), conn, params=exp_list)
df_IEJ['IEJ_id'] = df_IEJ['gene'] + "_" + df_IEJ['id'].astype(str)
df_IEJ_pivot=df_IEJ.pivot(index="IEJ_id", columns="exp", values="read_count")
df_IEJ_pivot=df_IEJ_pivot.fillna(0)

#make IEJ bed track for given exp
df_IEJbed = pd.read_sql("SELECT DISTINCT i.id,i.chr,x.start,x.end,i.strand,i.category,i.iso_exons,i.junctions, x.ex_sizes, x.ex_starts FROM isoform i INNER JOIN isoform_ends x on i.id = x.isoform_id WHERE i.IEJ='TRUE' AND x.exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
df_IEJbed['name'] = df_IEJbed['id'].map(str)+'_'+df_IEJbed['start'].map(str)+'_'+df_IEJbed['end'].map(str)
df_IEJbed['score'] = 60
df_IEJbed['RGB'] = df_IEJbed['category'].apply(lambda x: "29,47,95" if x == "full-splice_match" else '131,144,250' if x=="incomplete-splice_match" else '110,175,70' if x=="novel_in_catalog" else '250,199,72' if x=="novel_not_in_catalog" else '0,188,192' if x=="antisense" else '250,136,0' if x=="intergenic" else '172,0,0' if x=="genic" else '0,150,255' if x=='genic_intron' else '207,51,172' if x=="fusion" else '255,0,0')
df_IEJbed = df_IEJbed[['chr','start','end','name','score','strand','start','end','RGB','iso_exons','ex_sizes','ex_starts']]
df_IEJbed.to_csv("/mnt/data/chunlab/Christine/isoSeQL_TEST/ALS_2714_C030/20211011_IEJ.bed", sep='\t', index=False, header=False)
with open("20211011_IEJ.bed", "r+") as f: s=f.read(); f.seek(0); f.write("track visibility=2 itemRgb=On name='ALS2714_C030_IEJs'\n" + s)


#11/05/2021
#get list of isoforms that have many coordinates for the ends and how many ends
conn=sqlite3.connect("/mnt/data/chunlab/Christine/isoSeQL_TEST/2021nov_monoExonFix/20211105_2ALSsamp.db")
c=conn.cursor()
exp_list=[1]
diff_ends = pd.read_sql("SELECT isoform_id,COUNT(isoform_id) FROM isoform_ends WHERE exp IN (%s) GROUP BY isoform_id" % ",".join('?' for i in exp_list), conn, params=exp_list)
diff_ends=diff_ends.sort_values(by="COUNT(isoform_id)", ascending=False)


#by gene which isoforms have the most variable ends
gene_endsCounts = pd.read_sql("SELECT i.gene, x.isoform_id, COUNT(x.isoform_id) FROM isoform_ends x INNER JOIN isoform i on i.id=x.isoform_id WHERE x.exp IN (%s) GROUP BY x.isoform_id" % ','.join('?' for i in exp_list), conn, params=exp_list)
gene_endsCounts=gene_endsCounts.sort_values(by="COUNT(x.isoform_id)", ascending=False)

#which gene has the most (this isn't quite right b/c some genes repeated in isoform_ends)
gene_isCounts = pd.read_sql("SELECT i.gene,COUNT(i.gene) FROM isoform i INNER JOIN isoform_ends x on i.id=x.isoform_id WHERE x.exp IN (%s) GROUP BY gene" % ','.join('?' for i in exp_list), conn, params=exp_list)
gene_isCounts=gene_isCounts.sort_values(by="COUNT(i.gene)", ascending=False)

#get table of sample + exp info
info_table = pd.read_sql("SELECT e.id,s.id,e.seq_date,s.patient,s.tissue,s.disease,s.RIN,s.age FROM exp e INNER JOIN sampleData s on e.patient_id=s.id GROUP BY e.id",conn)
info_table.columns.values[0]="exp_id"
info_table.columns.values[1]="sample_id"
info_table.to_csv("/mnt/data/chunlab/Christine/isoSeQL_TEST/2021nov_monoExonFix/20211108_infoTable.txt", sep='\t', index=False, header=True)


#get table of experiment info
exp_table=pd.read_sql("SELECT * FROM exp e INNER JOIN sampleData s on e.patient_id=s.id GROUP BY e.id", conn)
exp_table.columns.values[0]="exp_id"
exp_table.columns.values[1]="sample_id"
exp_table = exp_table.drop(['id'], axis=1)


#11/12/2021
#counts matrix for tappAS - can make for just isoforms with same junctions regardless of ends and for isoforms with variable ends accounted for, input is explist
conn=sqlite3.connect("/mnt/data/chunlab/Christine/isoSeQL_TEST/2021nov_monoExonFix/20211105_2ALSsamp.db")
c=conn.cursor()
exp_list=[1,2]
commonJxn_counts = pd.read_sql("SELECT isoform_id, exp, read_count FROM counts WHERE exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
pivot=commonJxn_counts.pivot(index="isoform_id", columns="exp", values="read_count")
pivot=pivot.fillna(0)
#variable ends
varEnds_counts = pd.read_sql("SELECT isoform_id, start, end, exp, read_count FROM isoform_ends WHERE exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
varEnds_counts['isoform'] = varEnds_counts['isoform_id'].map(str)+'_'+varEnds_counts['start'].map(str)+'_'+varEnds_counts['end'].map(str)
pivot = varEnds_counts.pivot(index="isoform", columns="exp", values="read_count")
pivot=pivot.fillna(0)

#how to generate a gff3 file compatible with tappAS?
# PB.74447.1      tappAS  transcript      1       1119    .       -       .       ID=novel; primary_class=fusion; PosType=T
# PB.74447.1      tappAS  gene    1       1119    .       -       .       ID=CTD-2631K10.1_CTC-347C20.1; Name=CTD-2631K10.1_CTC-347C20.1; Desc=CTD-2631K10.1_CTC-347C20.1; PosType=T
# PB.74447.1      tappAS  CDS     .       .       .       -       .       ID=Protein_PB.74447.1; Name=Protein_PB.74447.1; Desc=Protein_PB.74447.1; PosType=T
# PB.74447.1      tappAS  genomic 1       1       .       -       .       Chr=chr5; PosType=G
# PB.74447.1      tappAS  splice_junction 72761497        72761614        .       -       .       ID=junction_1_canonical; Chr=chr5; PosType=G
# PB.74447.1      tappAS  splice_junction 72761671        72762609        .       -       .       ID=junction_2_canonical; Chr=chr5; PosType=G
# PB.74447.1      tappAS  splice_junction 72762925        72794398        .       -       .       ID=junction_3_canonical; Chr=chr5; PosType=G
# PB.74447.1      tappAS  splice_junction 72794485        72817249        .       -       .       ID=junction_4_canonical; Chr=chr5; PosType=G
#^^example of tappAS gff3

# gff3=pd.read_sql("SELECT DISTINCT i.id, i.chr, i.strand, i.gene, i.junctions, i.category, x.start, x.end, x.ex_sizes FROM isoform i INNER JOIN isoform_ends x on i.id=x.isoform_id WHERE x.exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
import numpy as np

gff3TX=pd.read_sql("SELECT DISTINCT i.id, i.chr, i.strand, i.gene, i.junctions, i.category, x.start, x.end, x.ex_sizes, t.tx FROM isoform i INNER JOIN isoform_ends x on i.id=x.isoform_id LEFT OUTER JOIN txID t ON i.id=t.isoform_id WHERE x.exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
gff3TX['isoform'] = gff3TX['id'].map(str)+'_'+gff3TX['start'].map(str)+'_'+gff3TX['end'].map(str)
gff3TX['tx'].replace(np.nan, "novel", inplace=True)
outFile = open("/mnt/data/chunlab/Christine/isoSeQL_TEST/2021nov_monoExonFix/tappAS_gff3.gff3", "w+")
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


#gene counts matrix
commonJxn_counts = pd.read_sql("SELECT isoform_id, exp, read_count FROM counts WHERE exp IN (%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)

gene_counts = pd.read_sql("SELECT c.isoform_id, c.exp, c.read_count, i.gene FROM counts c INNER JOIN isoform i on i.id=c.isoform_id WHERE exp IN(%s)" % ','.join('?' for i in exp_list), conn, params=exp_list)
test = gene_counts.groupby(['exp', 'gene'])['read_count'].sum().reset_index()


#tappASgff doesn't work b/c changed ends_id naming scheme
gff3TX=pd.read_sql("SELECT DISTINCT x.id, i.chr, i.strand, i.gene, i.junctions, i.category, x.start, x.end, x.ex_sizes, t.tx FROM isoform i LEFT OUTER JOIN txID t on i.id=t.isoform_id INNER JOIN isoform_ends x on x.isoform_id==i.id WHERE x.id IN (SELECT e.ends_id FROM ends_counts e WHERE e.exp IN (%s)) " % ','.join('?' for i in exp_list), conn, params=exp_list)

endsInfo=pd.read_sql("SELECT id, isoform_id,start, end, ex_sizes FROM isoform_ends WHERE id IN (SELECT ends_id FROM ends_counts WHERE exp IN (%s))" % ','.join('?' for i in exp_list), conn, params=exp_list)
isoInfo= pd.read_sql("SELECT i.id, i.chr, i.strand, i.gene, i.junctions, i.category, t.tx FROM isoform i INNER JOIN txID t on i.id=t.isoform_id WHERE i.id IN (SELECT isoform_id FROM counts WHERE exp IN (%s))" % ','.join('?' for i in exp_list), conn, params=exp_list)
isoInfo.rename(columns={'id':'isoform_id'}, inplace=True)
gff3TX=isoInfo.merge(endsInfo, on=['isoform_id'])

df_IEJ=pd.read_sql("SELECT i.id, i.gene, c.exp, c.read_count FROM isoform i INNER JOIN counts c on c.isoform_id = i.id WHERE i.IEJ = 'TRUE' AND c.exp IN (%s) " % ','.join('?' for i in exp_list), conn, params=exp_list)

df_IEJ=pd.read_sql("SELECT x.id, i.gene, c.exp, c.read_count FROM isoform_ends x INNER JOIN ends_counts c on x.id=c.ends_id INNER JOIN isoform i on x.isoform_id=i.id WHERE i.id IN (SELECT id FROM isoform WHERE IEJ='TRUE') AND c.exp IN (%s) " % ','.join("?" for i in exp_list), conn, params=exp_list)