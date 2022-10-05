#!/usr/bin/env Rscript

library(optparse)
library(UpSetR)

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
pdf(opt$outfile, useDingbats=FALSE, height=5, width=11)
print(upset(data, order.by="freq", nsets=num, nintersects=top, empty.intersections="on", mainbar.y.label="Isoforms in Common", sets.x.label="Isoforms per Group"), newpage=FALSE)
dev.off()