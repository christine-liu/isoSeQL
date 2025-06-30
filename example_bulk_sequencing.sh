#Example bulk sequencing analysis
#Author: Christine Liu

##Workman ONT samples (dRNA and cDNA)

SAMPLES=(NA12878-cDNA_All_Guppy_4.2.2 NA12878-DirectRNA_All_Guppy_4.2.2)

DATA_DIR=/path/to/data/
OUT_DIR=/path/to/output/
REF=/path/to/gencode48/GRCh38.primary_assembly.genome.fa


for SAMP in ${SAMPLES[@]}
do
	pbmm2 align --preset ISOSEQ -G 2500000 --sort ${REF} ${DATA_DIR}/${SAMP}.fastq ${OUT_DIR}/${SAMP}/${SAMP}_mapped_sorted.bam #map
	samtools view -b -q 10 -F 2304 -@ 16 -o ${OUT_DIR}/${SAMP}/${SAMP}_mapped_sorted_filtered.bam ${OUT_DIR}/${SAMP}/${SAMP}_mapped_sorted.bam #sort
	isoseq3 collapse --do-not-collapse-extra-5exons ${OUT_DIR}/${SAMP}/${SAMP}_mapped_sorted_filtered.bam  ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln.gff --min-aln-coverage 0.85 --min-aln-identity 0.8 #adjusted parameters to account for lower accuracy, collapse isoforms
done

#SQANTI3 annotation and filtering
ANNOT=/path/to/gencode48/gencode.v48.primary_assembly.annotation.gtf
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/sequence
for SAMP in ${SAMPLES[@]}
do
	cd ${OUT_DIR}/${SAMP}
	python /path/to/SQANTI3/sqanti3_qc.py ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln.gff ${ANNOT} ${REF} --genename --report skip --fl_count ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln.flnc_count.txt --skipORF
	python /path/to/SQANTI3/sqanti3_RulesFilter.py ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln_classification.txt ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln_corrected.fasta ${OUT_DIR}/${SAMP}/${SAMP}_mapped_filtered_collapsed_85cov80aln.gff --report skip
done

#add to isoSeQL database
DB_DIR=/path/to/isoSeQL_GM12878/
python /path/to/isoSeQL/isoSeQL_run.py --classif ${OUT_DIR}/NA12878-cDNA_All_Guppy_4.2.2/NA12878-cDNA_All_Guppy_4.2.2_mapped_filtered_collapsed_85cov80aln_classification.filtered_lite_classification.txt --genePred ${OUT_DIR}/NA12878-cDNA_All_Guppy_4.2.2/NA12878-cDNA_All_Guppy_4.2.2_mapped_filtered_collapsed_85cov80aln_corrected.genePred --sampleConfig ${DB_DIR}/samp_config/ONTcDNA_samp.config --expConfig ${DB_DIR}/exp_config/ONTcDNA_exp.config --db ${DB_DIR}/20250617_compare.db
python /path/to/isoSeQL/isoSeQL_run.py --classif ${OUT_DIR}/NA12878-DirectRNA_All_Guppy_4.2.2/NA12878-DirectRNA_All_Guppy_4.2.2_mapped_filtered_collapsed_85cov80aln_classification.filtered_lite_classification.txt --genePred ${OUT_DIR}/NA12878-DirectRNA_All_Guppy_4.2.2/NA12878-DirectRNA_All_Guppy_4.2.2_mapped_filtered_collapsed_85cov80aln_corrected.genePred --sampleConfig ${DB_DIR}/samp_config/ONTdRNA_samp.config --expConfig ${DB_DIR}/exp_config/ONTdRNA_exp.config --db ${DB_DIR}/20250617_compare.db

##GM12878_kinnex analysis

DATA_DIR=/path/to/data/
OUT_DIR=/path/to/output/

PRIMERS=/path/to/pb_barcodes/IsoSeq_v2_primers.fa
REF=/path/to/gencode48/GRCh38.primary_assembly.genome.fa
ANNOT=/path/to/gencode48/gencode.v48.primary_assembly.annotation.gtf
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/sequence

SAMP=GM12878
java -jar ~/picard.jar RevertSam -I ${DATA_DIR}/${SAMP}.kinnex.bam -O ${DATA_DIR}/${SAMP}.kinnex_unmapped.bam #needed unmapped bam to clip primer sequences
lima ${DATA_DIR}/${SAMP}.kinnex_umapped.bam ${PRIMERS} ${OUT_DIR}/${SAMP}_HIFI_demux.bam --isoseq --dump-clips
isoseq3 refine ${OUT_DIR}/${SAMP}_HIFI_demux.IsoSeqX_5p--IsoSeqX_3p.bam ${PRIMERS} ${OUT_DIR}/${SAMP}_HIFI_demux_flnc.bam  --require-polya 
isoseq cluster ${OUT_DIR}/${SAMP}_HIFI_demux_flnc.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster.bam --verbose --use-qvs
pbmm2 align --preset ISOSEQ -G 2500000 --sort ${REF} ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_mapped.bam
isoseq3 collapse --do-not-collapse-extra-5exons ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_mapped.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff
cd ${OUT_DIR}
python /path/to/SQANTI3/sqanti3_qc.py ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff ${ANNOT} ${REF} --genename --fl_count ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.flnc_count.txt --skipORF --report skip
python /path/to/SQANTI3/sqanti3_RulesFilter.py ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed_classification.txt ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed_corrected.fasta ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff --report skip
conda deactivate


source activate isoSeQL
DB_DIR=/path/to/isoSeQL_GM12878/
python /path/to/isoSeQL/isoSeQL_run.py --classif ${OUT_DIR}/GM12878_HIFI_demux_flnc_cluster_collapsed_classification.filtered_lite_classification.txt --genePred ${OUT_DIR}/GM12878_HIFI_demux_flnc_cluster_collapsed_corrected.genePred --sampleConfig ${DB_DIR}/samp_config/Kinnex_samp.config --expConfig ${DB_DIR}/exp_config/Kinnex_exp.config --db ${DB_DIR}/20250617_compare.db



##GM12878 encode analysis
DATA_DIR=/path/to/data/
OUT_DIR=/path/to/output/

PRIMERS=/path/to/pb_barcodes/isoSeqNEBprimers.fa
REF=/path/to/gencode48/GRCh38.primary_assembly.genome.fa
ANNOT=/path/to/gencode48/gencode.v48.primary_assembly.annotation.gtf
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/sequence

SAMP=ENCFF913VJX

ccs ${DATA_DIR}/ENCFF913VJX_subreads.bam ${OUT_DIR}/${SAMP}_HIFI.bam --report-file ${OUT_DIR}/${SAMP}_HIFI_report.txt
lima ${OUT_DIR}/${SAMP}_HIFI.bam ${PRIMERS} ${OUT_DIR}/${SAMP}_HIFI_demux.bam --isoseq --dump-clips
isoseq3 refine ${OUT_DIR}/${SAMP}_HIFI_demux.Clontech_5p--NEB_Clontech_3p.bam ${PRIMERS} ${OUT_DIR}/${SAMP}_HIFI_demux_flnc.bam  --require-polya 
isoseq cluster ${OUT_DIR}/${SAMP}_HIFI_demux_flnc.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster.bam --verbose --use-qvs
pbmm2 align --preset ISOSEQ -G 2500000 --sort ${REF} ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_mapped.bam
isoseq3 collapse --do-not-collapse-extra-5exons ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_mapped.bam ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff
cd ${OUT_DIR}
python /path/to/SQANTI3/sqanti3_qc.py ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff ${ANNOT} ${REF} --genename --fl_count ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.flnc_count.txt --skipORF --report skip
python /path/to/SQANTI3/sqanti3_RulesFilter.py ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed_classification.txt ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed_corrected.fasta ${OUT_DIR}/${SAMP}_HIFI_demux_flnc_cluster_collapsed.gff --report skip


DB_DIR=/path/to/isoSeQL_GM12878/
python /path/to/isoSeQL/isoSeQL_run.py --classif ${OUT_DIR}/ENCFF913VJX_HIFI_demux_flnc_cluster_collapsed_classification.filtered_lite_classification.txt --genePred ${OUT_DIR}/ENCFF913VJX_HIFI_demux_flnc_cluster_collapsed_corrected.genePred --sampleConfig ${DB_DIR}/samp_config/Encode_samp.config --expConfig ${DB_DIR}/exp_config/Encode_exp.config --db ${DB_DIR}/20250617_compare.db
