#!/bin/sh



### coverage per regions:


module spider mosdepth/0.3.1


for i in *.bam; do
mosdepth --by mono.5bp.final.bed \
--no-per-base --thresholds 3  "${i/%.bam/.mosdepth_3x.mono.5bp}"  "$i" 
done


rm *.mosdepth_3x.mono.5bp.regions.bed.gz
rm *.mosdepth_3x.mono.5bp.regions.bed.gz.csi
rm *.mosdepth_3x.mono.5bp.mosdepth.region.dist.txt
rm *.mosdepth_3x.mono.5bp.mosdepth.global.dist.txt

gunzip *.gz

for i in *.mosdepth_3x.thresholds.bed; do
cat "$i" | awk '{ SUM += $5} END { print SUM }' > "${i/%.thresholds.bed/.thresholds.txt}";
done



cut -f -3 *.mosdepth_3x.thresholds.txt > mev_CeNDR_mosdepth_3x.thresholds.txt

cut -f -3 *.mosdepth_3x.mono.5bp.thresholds.txt > mev_CeNDR_mosdepth_3x.mono.5bp.thresholds.txt

ls *.mosdepth_3x.thresholds.txt > list.txt


### 3mer_counts_c_elegans.WS263

module load jellyfish/2.3.0

jellyfish count -m 3 -s 101M -t 8 c_elegans.PRJNA13758.WS263.genomic.fa -o c_elegans.WS263.genomic.3mer.jf

jellyfish histo c_elegans.WS263.genomic.3mer.jf
 
jellyfish dump -c c_elegans.WS263.genomic.3mer.jf > 3mer_counts_c_elegans.WS263.genomic.txt



### 3bp_motif analysis:

module load bcftools
module load bedtools


# 0- extract SNPs:

for i in *.vcf; do
java -Xmx8g -jar /apps/gatk/4.1.4.0/gatk-package-4.1.4.0-local.jar SelectVariants \
    -R /ufrc/baer/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa \
    -V "$i" \
    --select-type-to-include SNP \
    -O "${i/%.vcf/.snp.vcf}"
done


# 1- extract list of positions:

for i in *.snp.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done


for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=0' "$i" > "${i/%.step2.bed/.step3.bed}"
done


for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done


# 2-  extract the sequence contex based on the bed file:

for i in *.final.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final.bed/.3bp_last.txt}"
done

# 3- counting the # of 3-bp motif:

for i in *.3bp_last.txt; do
    cut -f4 "$i" |  tr '[:space:]' '[\n*]' | grep -v "^\s*$" | sort | uniq -c | sort -bnr > "${i/%.3bp_last.txt/.3bp_results.txt}"
done



# 4- counting the number of each 3bp-motif separately:


for i in *.3bp_last.txt; do
    grep -o -i AAA "$i" | wc -l > "${i/%.3bp_last.txt/.3bp_AAA.rm.txt}"
done

cut -f1 *.3bp_AAA.rm.txt > AAA.ff.txt




