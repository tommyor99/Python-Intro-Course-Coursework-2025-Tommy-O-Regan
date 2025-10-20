grep -c ">" SzFaldDH_PFAM_blast_hits_filtered.fasta 
wc SzFaldDH_PFAM_blast_hits_filtered.fasta
head -400 SzFaldDH_PFAM_blast_hits_filtered.fasta > small_SzFaldDH_blast.fa
vi small_SzFaldDH_blast.fa
./aclust -p small_SzFaldDH_blast -jdis small_SzFaldDH_blast.fa 
./aclust -h | more
./bclust -p small_SzFaldDH_blast -jdis small_SzFaldDH_blast.fa 
ls -lt
ls -ltr
more small_SzFaldDH_blast.dbin.dat
bclust -h
./bclust -h
./bclust -p SzFaldDH_PFAM_blast_hits_filtered.fasta -jdis SzFaldDH_PFAM_blast_hits_filtered.fasta
history
exit
sudo apt install -y ncbi-blast+ mafft cd-hit seqkit python3-biopython build-essential git
git clone https://github.com/scapella/trimal.git
cd trimal/source
make
sudo cp trimal /usr/local/bin/
cd ../..  
blastp -version
mafft --version
cd-hit -h | head -n 2
seqkit version
trimal --help | head -n 3
trimal -h | head -n 20
efetch -db protein -id WP_010957160.1 -format fasta > seed.fasta
sudo apt update
sudo apt install ncbi-entrez-direct -y
efetch -db protein -id WP_010957160.1 -format fasta
sudo apt install hmmer
hmmsearch --tblout hits.tbl TIGRFAMs_15.0_HMMs/TIGR02819.hmm seed.fasta
efetch -db protein -id WP_043728606.1 -format fasta > seed.fasta
which efetch
cat seed.fasta
ls TIGRFAMs_15.0_HMMs
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz
tar -xvzf TIGRFAMs_15.0_HMM.tar.gz
ls
cd seed.fasta
ftp://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.tar.gz
ls
blastp -query seed.fasta -db nr -remote -max_target_seqs 500 -outfmt 6 -out blast_results.txt
ls
cut -f2 blast_results.txt | sort -u > hit_ids.txt
cat hit_ids.txt | epost -db protein | efetch -format fasta > family_raw.fasta
cad blast_results.txt 
cat blast_results.txt 
blastp -query seed.fasta -db nr -remote -max_target_seqs 500 -outfmt "6 sacc" -out blast_accs.txt
ls
sort -u blast_accs.txt > hit_ids.txt
cat hit_ids.txt | epost -db protein | efetch -format fasta > family_raw.fasta
cat hit_ids.txt 
ls
nano blast_accs.txt 
nano blast_results.txt 
git clone git@github.com:GarryGippert/Aclust.git
cd Aclust
cd src; make
ls
cd ..
ls
aclust -s ../dat/BLOSUM62.dat SzFaldDH_PFAM_blast_hits_filtered.fasta
