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


######### Making alignments that are just same patients
# Extract from processed alignment
for i in {1..68}; do seqkit grep -r -p "^${i}_" -p ON563414.3.masked /home/taouk/Monkeypox/alignments/consensus_pass_mafft_nogaps_masked_nodupref_0.1N.aln > ${i}.aln; done
ls *.aln > input.txt
for i in $(cat input.txt); do snipit ${i} -o ${i} -r ON563414.3.masked -f pdf; done


######### Global sequences making alignment from scratch
xz --decompress sequences.fasta.xz 
seqtk comp sequences.fasta > check.tsv ##check
dos2unix include.txt
seqkit grep -f include.txt sequences.fasta -o sequences_subset.fasta
cat /home/taouk/Monkeypox/alignments/consensus_pass.fasta sequences_subset.fasta > sequences_subset_aus.fasta
seqtk comp sequences_subset_aus.fasta > check2.tsv
#augur align --s sequences_subset_aus.fasta --o sequences_subset_aus.aln --nthreads 40 --reference-name ON563414.3 
cd Monkeypox/nextclade
conda activate nextclade
nextclade dataset get --name hMPXV --output-dir hmpxv
nextclade run -D hmpxv/ /home/taouk/Monkeypox/global/sequences_subset_aus.fasta --output-all local
cp /home/taouk/Monkeypox/nextclade/global/nextclade.aligned.fasta /home/taouk/Monkeypox/global/nextclade.aligned.fasta
cd /home/taouk/Monkeypox/global
conda activate Base
trimal -in nextclade.aligned.fasta -out nextclade.aligned_nogaps.fasta -gt 1
seqkit grep -p "ON563414.3.masked" nextclade.aligned_nogaps.fasta -o ON563414.3.masked.fasta
dos2unix mask.txt
for i in $(cat include.txt); do sed "s/chr/${i}/g" mask.txt >> mask_all.bed; done
bedtools maskfasta -fi nextclade.aligned_nogaps.fasta -fo nextclade.aligned_nogaps_masked.fasta -bed mask_all.bed
/home/gtaiaroa/Capture/Paper/Align/Scripts/goalign-master/goalign clean sites --char=N -c 0.1 --ignore-case -i nextclade.aligned_nogaps_masked.fasta -o nextclade.aligned_nogaps_masked_N0.1.fasta
seqkit grep -vp "ON563414.3.masked" nextclade.aligned_nogaps_masked.fasta -o nextclade.aligned_nogaps_masked_nodupref.fasta
seqkit grep -vp "ON563414.3.masked" nextclade.aligned_nogaps_masked_N0.1.fasta -o nextclade.aligned_nogaps_masked_N0.1_nodupref.fasta
seqkit grep -vp "ON563414.3.masked" nextclade.aligned_nogaps.fasta -o nextclade.aligned_nogaps_nodupref.fasta
#Trees
cd Tree
iqtree -s nextclade.aligned_nogaps_masked_nodupref.fasta -B 1000 -nt 40 -m HKY+F+I
cd ../Tree2
iqtree -s nextclade.aligned_nogaps_masked_nodupref.fasta -B 1000 -nt 30 --polytomy -m HKY+F+I
cd ../Tree3
cp ../nextclade.aligned_nogaps_nodupref.fasta nextclade.aligned_nogaps_nodupref.fasta
iqtree -s nextclade.aligned_nogaps_nodupref.fasta -B 1000 -nt 30 --polytomy -m HKY+F+I


######### Nextclade
cd Monkeypox/nextclade
seqkit grep -vrp "ON563414.3" /home/taouk/Monkeypox/alignments/consensus_pass.fasta > local.fasta
conda activate nextclade

#|hMPXV                              │ NC_063383.1 (*)                   │ 2022-09-27T12:00:00Z (*) │ name=hMPXV                    │ Add lineages A.2.2 and B.1.10-12  │
#│ 'Human Monkeypox (hMPXV)'         │ 'MPXV-M5312_HM12_Rivers'          │                          │ reference=NC_063383.1 (*)     │                                   │
#│                                   │                                   │                          │ tag=2022-09-27T12:00:00Z (*)  │                      

nextclade dataset get --name hMPXV --output-dir hmpxv
nextclade run -D hmpxv/ local.fasta --output-all local





