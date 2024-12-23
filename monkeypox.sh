#!/bin/bash

##### Mona Taouk ######

#######################
########### ###########
###### MONKEYPOX ######
########### ###########
#######################


#########Testing thresholds for Twist data
cat /home/taouk/Monkeypox/ont/ont_consensus/32_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/0/32_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/3/32_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/60/32_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/75/32_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/90/32_A.consensus.fasta > /home/taouk/Monkeypox/Compare_thresholds/32_A.fasta
cat /home/taouk/Monkeypox/ont/ont_consensus/60_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/0/60_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/3/60_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/60/60_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/75/60_A.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/90/60_A.consensus.fasta > /home/taouk/Monkeypox/Compare_thresholds/60_A.fasta
cat /home/taouk/Monkeypox/ont/ont_consensus/60_D.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/0/60_D.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/3/60_D.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/60/60_D.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/75/60_D.consensus.fasta /home/taouk/Monkeypox/twist/twist_consensus/90/60_D.consensus.fasta > /home/taouk/Monkeypox/Compare_thresholds/60_D.fasta
mafft --thread 40 /home/taouk/Monkeypox/Compare_thresholds/32_A.fasta > /home/taouk/Monkeypox/Compare_thresholds/32_A.aln
mafft --thread 40 /home/taouk/Monkeypox/Compare_thresholds/60_A.fasta > /home/taouk/Monkeypox/Compare_thresholds/60_A.aln
mafft --thread 40 /home/taouk/Monkeypox/Compare_thresholds/60_D.fasta > /home/taouk/Monkeypox/Compare_thresholds/60_D.aln


######### Masking reference
bedtools maskfasta -fi mpx_us_22.fasta -fo mpx_us_22_masked.fasta -bed MPXV_mask.bed


######### multifasta of illumina and ONT pass
for i in $(cat twist_90_pass.txt); do cat /home/taouk/Monkeypox/twist/twist_consensus/75/${i}.consensus.fasta >> consensus_pass_oldnames.fasta; done
for i in $(cat ont_85_pass.txt); do cat /home/taouk/Monkeypox/ont/ont_consensus/${i}.consensus.fasta >> consensus_pass_oldnames.fasta; done
#rename headers
seqkit replace -p '^(\S+)' -r '{kv}' -k rename.txt consensus_pass_oldnames.fasta > consensus_pass.fasta
#add refs
cat /home/taouk/Monkeypox/mpx_us_22.fasta /home/taouk/Monkeypox/mask_ref/mpx_us_22_masked.fasta >> consensus_pass.fasta
#check
seqtk comp consensus_pass.fasta > stats.tsv
rm consensus_pass_oldnames.fasta


######### Making alignments
mafft --thread 40 consensus_pass.fasta > consensus_pass_mafft.aln
# Remove all gaps from alignment (removes columns with any gaps)
trimal -in consensus_pass_mafft.aln -out consensus_pass_mafft_nogaps.aln -gt 1
# Make bedfile
for i in $(cat IDs.txt); do sed "s/chr/${i}/g" mask_nogaps.txt >> mask_nogaps_all.bed; done
# Masking alignment with new bedfile
bedtools maskfasta -fi consensus_pass_mafft_nogaps.aln -fo consensus_pass_mafft_nogaps_masked.aln -bed mask_nogaps_all.bed
# Remove that other reference
seqkit grep -vp ON563414.3 consensus_pass_mafft_nogaps_masked.aln > consensus_pass_mafft_nogaps_masked_nodupref.aln
# Removing Ns
/home/gtaiaroa/Capture/Paper/Align/Scripts/goalign-master/goalign clean sites --char=N -c 0 --ignore-case -i consensus_pass_mafft_nogaps_masked_nodupref.aln -o consensus_pass_mafft_nogaps_masked_nodupref_noNs.aln
/home/gtaiaroa/Capture/Paper/Align/Scripts/goalign-master/goalign clean sites --char=N -c 0.1 --ignore-case -i consensus_pass_mafft_nogaps_masked_nodupref.aln -o consensus_pass_mafft_nogaps_masked_nodupref_0.1N.aln
# ML Tree
iqtree -s consensus_pass_mafft_nogaps_masked_nodupref_0.1N.aln -B 1000 -nt 40 --polytomy
#snp distances
snp-dists consensus_pass_mafft_nogaps_masked_nodupref_0.1N.aln > snp_dist.tsv


######### Nextclade Local
cd Monkeypox/nextclade
seqkit grep -vrp "ON563414.3" /home/taouk/Monkeypox/alignments/consensus_pass.fasta > local.fasta
conda activate nextclade
nextclade dataset list

#|hMPXV                              │ NC_063383.1 (*)                   │ 2022-09-27T12:00:00Z (*) │ name=hMPXV                    │ Add lineages A.2.2 and B.1.10-12  │
#│ 'Human Monkeypox (hMPXV)'         │ 'MPXV-M5312_HM12_Rivers'          │                          │ reference=NC_063383.1 (*)     │                                   │
#│                                   │                                   │                          │ tag=2022-09-27T12:00:00Z (*)  │                      


######### Nextclade Global
conda activate Base
cd Monkeypox/nextclade
dos2unix include.txt
seqkit grep -f include.txt /home/taouk/Monkeypox/global/PREVIOUS/sequences.fasta > global.fasta
dos2unix include_new.txt
seqkit grep -f include_new.txt /home/gtaiaroa/Mona/Victoria/Set.fasta > new_SNP.fasta
cat local.fasta global.fasta new_SNP.fasta /home/taouk/Monkeypox/mpx_us_22.fasta > sequences_global_local.fasta
seqtk comp sequences_global_local.fasta > comp.tsv
conda activate nextclade
nextclade run -D hmpxv/ sequences_global_local.fasta --output-all global_new
conda activate Base
cd global_new
snp-dists nextclade.aligned.fasta > nextclade.aligned.distances.tsv


######### Global Trees
cd global
cp /home/taouk/Monkeypox/nextclade/global_new/nextclade.aligned.fasta nextclade.aligned.fasta
cd /home/taouk/Monkeypox/global
conda activate Base
trimal -in nextclade.aligned.fasta -out nextclade.aligned_nogaps.fasta -gt 1
/home/gtaiaroa/Capture/Paper/Align/Scripts/goalign-master/goalign clean sites --char=N -c 0 --ignore-case -i nextclade.aligned_nogaps.fasta -o nextclade.aligned_nogaps_noNs.fasta
#Trees
cd global_newSNP
cp /home/taouk/Monkeypox/global/nextclade.aligned.fasta nextclade.aligned.fasta
conda activate Base
iqtree -s nextclade.aligned.fasta -B 1000 -nt 40 -m HKY+F+I --polytomy


# Global tree with position missing
cd Monkeypox/D1604N/new_alignment
iqtree -s Deleted_new.fasta -B 1000 -nt 40 --polytomy -m HKY+F+I


######### Alignments per patient using Rivers reference and masking for repeats (attempt)
conda activate nextclade
cd /home/taouk/Monkeypox/patient_alignment/rivers_usa
nextclade run -D /home/taouk/Monkeypox/nextclade/hmpxv sequences.fasta --output-all local
# seperate
conda activate Base
for i in {1..68}; do seqkit grep -r -p "^${i}_" -p ON563414.3.masked  /home/taouk/Monkeypox/patient_alignment/rivers_usa/nextclade.aligned.fasta > individuals/${i}.aln; done
cd individuals
ls *.aln > input.txt
for i in $(cat input.txt); do snipit ${i} -o ${i} -r ON563414.3.masked -f pdf --height 3 --width 10; done
######### Alignments removing the reference
for i in {1..68}; do seqkit grep -pv 


######### SUmmary of Snipit for whole aligment
snipit nextclade.aligned_edit.fasta -o all -r ON563414.3.masked -f pdf
snp-sites -v -o test nextclade.aligned_edit.fasta


####### VCFs
/home/taouk/Monkeypox/ont/ont_consensus/66_D.pass.vcf
/home/taouk/Monkeypox/ont/ont_consensus/66_C.pass.vcf.gz
/home/taouk/Monkeypox/ont/ont_consensus/15_B.pass.vcf.gz
/home/taouk/Monkeypox/ont/ont_consensus/15_A.pass.vcf.gz





