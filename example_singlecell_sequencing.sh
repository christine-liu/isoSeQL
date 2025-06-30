#Example single-cell sequencing

#start with tagged, refined, corrected.sorted.dedup.bam that has been group dedup'd


REF=/path/to/gencode44/GRCh38.p14.genome.fa
ANNOT=/path/to/gencode44/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
TENK=/data/cliu/PacBio_Kinnex/SQ2_PBMC_10k/ #10k sample data
FIVEK=/data/cliu/PacBio_Kinnex/SQ2_PBMC_5k/ #5k sample data


##download fastas
samtools fasta ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam > ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta
samtools fasta ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam > ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta
samtools view -h -o ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.sam ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.bam
samtools view -h -o ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.sam ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.bam

#collapse isoforms
collapse_isoforms_by_sam.py --input ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta -s ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.sam -c 0.99 -i 0.95 --gen_mol_count -o ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake --cpus 8
collapse_isoforms_by_sam.py --input ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta -s ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_mapped.sam -c 0.99 -i 0.95 --gen_mol_count -o ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake --cpus 8

#SQANTI3
export PYTHONPATH=$PYTHONPATH:/path/to/cDNA_Cupcake/sequence

cd ${TENK}
python /path/to/SQANTI3/sqanti3_qc.py scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.gff ${ANNOT} ${REF} --genename --fl_count scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.abundance.txt --skipORF --report skip
python /path/to/SQANTI3/sqanti3_RulesFilter.py scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_classification.txt scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_corrected.fasta scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_corrected.gtf --report skip

cd ${FIVEK}
python /path/to/SQANTI3/sqanti3_qc.py scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.gff ${ANNOT} ${REF} --genename --fl_count scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.abundance.txt --skipORF --report skip
python /path/to/SQANTI3/sqanti3_RulesFilter.py scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_classification.txt scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_corrected.fasta scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_corrected.gtf --report skip
conda deactivate


python /path/to/cDNA_Cupcake/singlecell/make_csv_for_dedup.py -i ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta -o ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv
gzip ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv
python /path/to/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.group.txt ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv.gz ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_classification.filtered_lite_classification.txt ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv
gzip ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv

python /path/to/cDNA_Cupcake/singlecell/make_csv_for_dedup.py -i ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta -o ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv
gzip ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv
python /path/to/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed.group.txt ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_info.csv.gz ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_cupcake.collapsed_classification.filtered_lite_classification.txt ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv
gzip ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv

#reverse complement barcodes
awk -F $',' '
     BEGIN {
         old="ATGC"
         new="TACG"
         for (i=1;i<=length(old);i++) {
             tr[substr(old,i,1)] = substr(new,i,1)
         }
     }
     {
         newVar=""
         for (i=1;i<=length($1);i++) {
             char = substr($1,i,1)
             newVar = (char in tr ? tr[char] : char) newVar
         }
         print newVar,$2
     }' OFS=',' ${FIVEK}/Seurat/SQ2_5k_celltypes >> ${FIVEK}/Seurat/SQ2_5k_celltypes_rc
python /path/to/cDNA_Cupcake/singlecell/link_molecule_to_celltype.py ${FIVEK}/Seurat/SQ2_5k_celltypes_rc ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv.gz ${FIVEK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated_celltype.csv

awk -F $',' '
     BEGIN {
         old="ATGC"
         new="TACG"
         for (i=1;i<=length(old);i++) {
             tr[substr(old,i,1)] = substr(new,i,1)
         }
     }
     {
         newVar=""
         for (i=1;i<=length($1);i++) {
             char = substr($1,i,1)
             newVar = (char in tr ? tr[char] : char) newVar
         }
         print newVar,$2
     }' OFS=',' ${TENK}/Seurat/SQ2_10k_celltypes >> ${TENK}/Seurat/SQ2_10k_celltypes_rc

python /path/to/cDNA_Cupcake/singlecell/link_molecule_to_celltype.py ${TENK}/Seurat/SQ2_10k_celltypes_rc ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated.csv.gz ${TENK}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup_annotated_celltype.csv



