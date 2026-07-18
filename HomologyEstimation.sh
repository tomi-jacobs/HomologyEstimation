sed -r -i 's/^>TRINITY/>ANLI2/' ANLI2.transdecoder.cds.fasta | head

cat *.cds.fasta > 39.fasta
makeblastdb -in 39.fasta -dbtype nucl -out 39.fasta
blastn -db 39.fasta -query 39.fasta -evalue 10 -num_threads 20 -max_target_seqs 1000 -out 39.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'
conda activate python2_env
python2 ~/data/YaYang/master/scripts/blast_to_mcl.py 39.rawblast 0.3
mcl 39.rawblast.hit-frac0.3.minusLogEvalue --abc -te 4 -tf 'gq(5)' -I 1.4 -o hit-frac0.3_I1.4_e5
mkdir fastaFilesMCL
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_mcl.py ../39.fasta ../hit-frac0.3_I1.4_e5 4 .

Round 1
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
parallel iqtree2 \-s {} \-nt 2 \-m GTR+G \-bb 1000 ::: *-cln
sed -i "s/\_D/\@D/g" *.treefile
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
&& runs two lines of code simultaneously 
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
• We are here right now with Joe
Renaming: sed -i "s/\_D/\@/g" *-cln
mkdir subtrees1
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py .. .mm 0.5 4 ../subtrees1/
Deviation for prelim for orthologs
mkdir Orthologs
python2 ~/data/YaYang/master/scripts/prune_paralogs_MI.py subtrees1/ .subtree 0.4  0.5 20 Orthologs/
ls Orthologs/ | wc -l
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py ../../../39.fasta .. .subtree ../fastaFromSubtrees1/

python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round3/subtrees3/fastafromFiles3   --r2 ~/data/fasta/prelim/round4/subtrees4/fastafromFiles4   --out_cluster 4cluster_membership_and_sequences.csv   --out_changed_only 4clusters_changed_only.csv   --out_taxon_counts 4taxon_counts_by_cluster.csv

python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 1cluster_membership_and_sequences.csv
End. Repeat.
 
Round 2
mkdir round 2
cp *.fa -r ~/data/fasta/prelim/round1/subtrees1/fastaFromSubtrees1/ ~/data/fasta/prelim/round2/
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round2/ .mm 0.5 4 ~/data/fasta/prelim/round2/substrees2/
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py ../../39.fasta ~/data/fasta/prelim/round2/subtrees2/ .subtree ~/data/fasta/prelim/round2/subtrees2/fastafromFiles2/

Round 3
cp ~/data/fasta/prelim/round2/subtrees2/fastafromFiles2/Change/*.fa ~/data/fasta/prelim/round3/
ls round3/ | wc –l (6553)
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
conda activate python2_env
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round3/ .mm 0.5 4 ~/data/fasta/prelim/round3/subtrees3/
sed -i "s/\_DN/\@DN/g" *.subtree
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round3/subtrees3/ .subtree ~/data/fasta/prelim/round3/subtrees3/fastafromFiles3/

Code for differences
python3 gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round1/subtrees1/fastafromFiles   --r2 ~/data/fasta/prelim/round2/subtrees2/fastafromFiles2   --out_cluster cluster_membership_and_sequences.csv   --out_changed_only clusters_changed_only.csv   --out_taxon_counts taxon_counts_by_cluster.csv

# 2) Run it (replace paths with your actual round-1 and round-2 folders):
python3 compare_rounds.py --r1 /path/to/round1_fastas --r2 /path/to/round2_fastas

python3 gene_membership_only.py \
  --r1 ~/data/fasta/prelim/round1/subtrees1/fastafromFiles \
  --r2 ~/data/fasta/prelim/round2/subtrees2/fastafromFiles2 \
  --out_cluster cluster_membership_changes.csv \
  --out_changed_only clusters_changed_only.csv \
  --out_taxon taxon_membership_summary.csv

python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round2/subtrees2/fastafromFiles2   --r2 ~/data/fasta/prelim/round3/subtrees3/fastafromFiles3   --out_cluster 3cluster_membership_and_sequences.csv   --out_changed_only 3clusters_changed_only.csv   --out_taxon_counts 3taxon_counts_by_cluster.csv
The above worked

python3 MyScript.py sorted_cluster_membership_and_sequences.csv
python3 Change_NoChange_Joe.py 3cluster_membership_and_sequences.csv


python3 ~/data/fasta/prelim/round1/subtrees1/fastafromFiles/r1_r2_diff_report_2.py
~/data/fasta/prelim/round3/subtrees3/fastafromFiles3

How many .fa in a directory? ls -1 *.fa 2>/dev/null | wc -l

Round 4
cp ~/data/fasta/prelim/round3/subtrees3/fastafromFiles3/Change/*.fa ~/data/fasta/prelim/round4/
ls round3/ | wc –l (717)
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round4/ .mm 0.5 4 ~/data/fasta/prelim/round4/subtrees4/
#python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round4/ .mm 0.5 4 ~/data/fasta/prelim/round4/
mkdir fastafromFiles4 
sed -i "s/\_DN/\@DN/g" *.subtree
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round4/subtrees4/ .subtree ~/data/fasta/prelim/round4/subtrees4/fastafromFiles4/
python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round3/subtrees3/fastafromFiles3   --r2 ~/data/fasta/prelim/round4/subtrees4/fastafromFiles4   --out_cluster 4cluster_membership_and_sequences.csv   --out_changed_only 4clusters_changed_only.csv   --out_taxon_counts 4taxon_counts_by_cluster.csv
python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 4cluster_membership_and_sequences.csv

Round 5
cp ~/data/fasta/prelim/round4/subtrees4/fastafromFiles4/Change/*.fa ~/data/fasta/prelim/round5/
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round5/ .mm 0.5 4 ~/data/fasta/prelim/round5/subtrees5/
sed -i "s/\_DN/\@DN/g" *.subtree
Mkdir fastafromFiles5
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round5/subtrees5/ .subtree ~/data/fasta/prelim/round5/subtrees5/fastafromFiles5/
python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round4/subtrees4/fastafromFiles4   --r2 ~/data/fasta/prelim/round5/subtrees5/fastafromFiles5   --out_cluster 5cluster_membership_and_sequences.csv   --out_changed_only 5clusters_changed_only.csv   --out_taxon_counts 5taxon_counts_by_cluster.csv
python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 5cluster_membership_and_sequences.csv

Round 6
cp ~/data/fasta/prelim/round5/subtrees5/fastafromFiles5/Change/*.fa ~/data/fasta/prelim/round6/
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round6/ .mm 0.5 4 ~/data/fasta/prelim/round6/subtrees6/
sed -i "s/\_DN/\@DN/g" *.subtree
Mkdir fastafromFiles6
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round6/subtrees6/ .subtree ~/data/fasta/prelim/round6/subtrees6/fastafromFiles6/
python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round5/subtrees5/fastafromFiles5   --r2 ~/data/fasta/prelim/round6/subtrees6/fastafromFiles6   --out_cluster 6cluster_membership_and_sequences.csv   --out_changed_only 6clusters_changed_only.csv   --out_taxon_counts 6taxon_counts_by_cluster.csv
python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 6cluster_membership_and_sequences.csv

Round 7
cp ~/data/fasta/prelim/round6/subtrees6/fastafromFiles6/Change/*.fa ~/data/fasta/prelim/round7/
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round7/ .mm 0.5 4 ~/data/fasta/prelim/round7/subtrees7/
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round7/subtrees7/ .subtree ~/data/fasta/prelim/round7/subtrees7/fastafromFiles7/
python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round6/subtrees6/fastafromFiles6   --r2 ~/data/fasta/prelim/round7/subtrees7/fastafromFiles7   --out_cluster 7cluster_membership_and_sequences.csv   --out_changed_only 7clusters_changed_only.csv   --out_taxon_counts 7taxon_counts_by_cluster.csv
python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 7cluster_membership_and_sequences.csv

Round 8
cp ~/data/fasta/prelim/round7/subtrees7/fastafromFiles7/Change/*.fa ~/data/fasta/prelim/round8/
taskset --cpu-list 0-23 parallel mafft \-\-auto \-\-maxiterate 1000 \-\-thread 2 {} \> {}\.aln ::: *.fa
for x in *.aln; do ~/data/phyx/src/pxclsq -s $x -p 0.1 -o $x-cln; done
taskset --cpu-list 0-20 parallel iqtree2 \-s {} \-nt 4 \-m GTR+G \-bb 1000 ::: *-cln
python2 ~/data/YaYang/master/scripts/trim_tips.py . .treefile 0.1 0.2
python2 ~/data/YaYang/master/scripts/mask_tips_by_taxonID_transcripts.py . . y
Mkdir subtrees8
python2 ~/data/YaYang/master/scripts/cut_long_internal_branches.py ~/data/fasta/prelim/round8/ .mm 0.5 4 ~/data/fasta/prelim/round8/subtrees8/
sed -i "s/\_DN/\@DN/g" *.subtree
python2 ~/data/YaYang/master/scripts/write_fasta_files_from_trees.py /home/tomi/data/fasta/prelim/39.fasta ~/data/fasta/prelim/round8/subtrees8/ .subtree ~/data/fasta/prelim/round8/subtrees8/fastafromFiles8/
python3 ~/data/projectOne/scripts/gene_membership_and_sequences.py   --r1 ~/data/fasta/prelim/round7/subtrees7/fastafromFiles7   --r2 ~/data/fasta/prelim/round8/subtrees8/fastafromFiles8   --out_cluster 8cluster_membership_and_sequences.csv   --out_changed_only 8clusters_changed_only.csv   --out_taxon_counts 8taxon_counts_by_cluster.csv
python3 ~/data/projectOne/scripts/Change_NoChange_Joe.py 8cluster_membership_and_sequences.csv
