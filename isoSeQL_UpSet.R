#!/usr/bin/env Rscript

library(optparse)
# library(UpSetR)
library(ComplexUpset)
library(ggplot2)

option_list=list(
	make_option(c('-i', '--input'), help='input binary matrix with experiments as columns'),
	make_option(c('-o', '--outfile'), help='path to output file'),
	make_option(c('-n', '--numSamples'), help='number of groups being compared'),
	make_option(c('-t', '--top'), help='show top x overlaps')
)
opt= parse_args(OptionParser(option_list=option_list))

data=read.csv(opt$input, header=TRUE, sep='\t')
num=opt$numSamples
top=opt$top
pdf(opt$outfile, useDingbats=FALSE, height=as.integer(num)/1.5+4, width=as.integer(top)/2)
# print(upset(data, order.by="freq", nsets=num, nintersects=top, empty.intersections="on", mainbar.y.label="Isoforms in Common", sets.x.label="Isoforms per Group"), newpage=FALSE)
groups=colnames(data)[3:(as.integer(num)+2)]
data[groups]=data[groups]==1
print(upset(data, groups, n_intersections=as.integer(top), base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background="black", on_bar="black"), mapping=aes(fill=factor(category, levels=c("fusion", "genic", "intergenic", "antisense", "novel_not_in_catalog", "novel_in_catalog", "incomplete-splice_match", "full-splice_match"))))+scale_fill_manual(values=c("fusion"="#CF33AC", "genic"="#AC0000", "intergenic"="#FA8800", "antisense"="#00BCC0", "novel_not_in_catalog"="#FAC748", "novel_in_catalog"="#6EAF46", "incomplete-splice_match"="#8390FA","full-splice_match"="#1D2F5F"))+labs(fill='Category')), width_ratio=0.25, set_sizes=(upset_set_size() + theme(axis.text.x=element_text(angle=90)))))
dev.off()