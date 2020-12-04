# My-genes-analysis-
 NW_001834246.1
Create a Fasta file of gene
ncbi-acc-download -F fasta -m protein XP_001628345.2 Now, perform a blast search using the query protein

Blast Query protein
blastp -db ~/data/blast/allprotein.fas -query XP_001628345.2.fa -outfmt 0 -max_hsps 1 > XP_001628345.blastp.typical.out
blastp -db ~/data/blast/allprotein.fas -query XP_001628345.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out XP_001628345.blastp.detail.out 

Filter blast
awk '{if ($6<0.00000000000001)print $1 }' XP_001628345.blastp.detail.out > XP_001628345.blastp.detail.filtered.out
 
Align Gene Family
seqkit grep --pattern-file XP_001628345.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_001628345.blastp.detail.filtered.fa
muscle -in XP_001628345.blastp.detail.filtered.fas -out XP_001628345.blastp.detail.filtered.aligned.fas
 t_coffee -other_pg seq_reformat -in XP_001628345.blastp.detail.filtered.aligned.fas -output simt_coffee -other_pg seq_reformat -in XP_001628345.blastp.detail.filtered.aligned.fas -output sim

Removed gaps t_coffee
t_coffee -other_pg seq_reformat -in XP_001628345.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out allhomologs.aligned.r50.fa

Constructing a Phylogenetic tree
iqtree -s XP_001628345.aligned.r50.fas -nt 2
gotree draw png -w 1000 -i XP_001628345.aligned.r50.fas.treefile -r -o XP_001628345.aligned.r50.fas.png

midpoint rooting
gotree reroot midpoint -i XP_001628345.aligned.r50.fas.treefile -o XP_001628345.aligned.r50.fas.midpoint.treefile
draw tree
nw_display -s XP_001628345.aligned.r50.fas.midpoint.treefile -w 1000 -b 'opacity:0' > XP_001628345.aligned.r50.fas.midpoint.svg
nw_reroot XP_001628345.aligned.r50.fas.treefile >XP_001628345.aligned.r50.fas.MendozaRoot.treefile
XP_001628345.aligned.r50.fas.MendozaRoot.treefile -w 1000 -b 'opacity:0' > XP_001628345.aligned.r50.fas.MendozaRoot.svg

cladogram
XP_001628345.aligned.r50.fas.MendozaRoot.treefile | nw_display -s -w 1000 > XP_001628345.aligned.r50.fas.MendozaRoot.cladogram.svg -

Bootstrap
iqtree -s XP_001628345.aligned.r50.fas -bb 1000 -nt 2 --prefix XP_001628345.r50.ufboot
nw_reroot XP_001628345.r50.ufboot.treefile Drosophila_melanogaster_Disks_large_1_DLG1_P31007 Nematostella_vectensis_disks_large_homolog_1_XP_001638123.2 Pocillopora_damicornis_A0A3M6TL78 Homo_sapiens_Disks_large_homolog_3_DLG3_Q92796 Homo_sapiens_Disks_large_homolog_4_DLG4_P78352 Homo_sapiens_Disks_large_homolog_1_DLG1_Q12959 Homo_sapiens_Disks_large_homolog_2_DLG2_Q15700 >XP_001628345.r50.ufboot.midpoint.treefile 
nw_display -s XP_001628345.r50.ufboot.midpoint.treefile -w 1000 -b 'opacity:0' > XP_001628345.r50.ufboot.midpoint.svg

Reconciliation
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar --help  
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b toybatch.txt --reconcile --speciestag prefix --savepng --treestats --events --homologtabletabs --phylogenomics 
sudo easy_install -U ete3 
python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g toygenetree.tre.reconciled --include.species
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b zobatch.txt --reconcile --speciestag prefix  --savepng --treestats --events  --phylogenomics 
python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g XP_001628345.r50.ufboot.midpoint.treefile --include.species
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b zobatch.txt --root --speciestag prefix  --savepng --treestats --events  --phylogenomics 

Rearrangement via Notung
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b mybatch.txt --rearrange --speciestag prefix --savepng --treestats --events --outputdir mygene --edgeweights name --threshold 90
python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g mygene/XP_001628345.r50.ufboot.midpoint.treefile.rearrange.0 --include.species

Re-arranged IQ-tree
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -g zoreconcileRearrange/XP_001630613.r50.ufboot.MendozaRoot.treefile.rearrange.0   -s species.tre --reconcile --speciestag prefix  --treeoutput newick --nolosses
gotree unroot -i XP_001628345.r50.ufboot.midpoint.treefile.rearrange.0.reconciled -o XP_001628345.r50.ufboot.unrooted.treefile.rearrange
cat XP_001628345.r50.ufboot.treefile XP_001628345.r50.ufboot.unrooted.treefile.rearrange > XP_001628345.r50.alternativetrees
iqtree -s XP_001628345.aligned.r50.fa -z XP_001628345.r50.alternativetrees -au -zb 10000 --prefix ZO1_altTrees -m LG+F+R5 -nt 2 -te XP_001628345.r50.ufboot.treefile

Domain Identification
iprscan5   --email Fatou.jobe@stonybrook.edu  --multifasta --useSeqId --sequence   XP_001628345.blastp.detail.filtered.ren.fas
cat *.tsv.txt > maguk.domains.all.tsv
grep Pfam maguk.domains.all.tsv >  maguk.domains.pfam.tsv
awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' maguk.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > maguk.domains.pfam.evol.tsv
 cut -f 1 maguk.domains.pfam.tsv | sort | uniq -c
iprscan5   --email your.email@stonybrook.edu  --multifasta --useSeqId --sequence   yourgenefamily.blastp.detail.filtered.ren.fas
cat *.tsv.txt > yourgenefamily.domains.all.tsv
grep Pfam yourgenefamily.domains.all.tsv >  yourgenefamily.domains.pfam.tsv
 awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' yourgenefamily.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > yourgenefamily.domains.pfam.evol.tsv

