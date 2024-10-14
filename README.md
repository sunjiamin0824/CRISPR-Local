CRISPR-Local
============

Please forward any question and suggestion about CRISPR-Local to: 
## sunjiamin0824@qq.com
The original email address 'sunjm0824@webmail.hzau.edu.cn' has been deactivated.

## A local single-guide RNA (sgRNA) design tool for non-reference plant genomes.
### Motivation:
For genome editing experiments, it is critical to design a reliable single-guide RNA (sgRNA). However, most existing tools are limited to reference genome and do not consider variations in genetic background. Taking this variation into account is especially important for plant genome editing studies, since the transformed lines are rarely the sequenced reference genomes and the presence of large diversity between them.

### Highlights:




* Designing sgRNAs suitable for non-reference varieties. In the plants for gene editing, the genome of the material and the reference genome are quite different. Leading to a fact that the sgRNAs designed according to the reference genome can be not well used for other materials. 
* Supporting various PAM. PAM motifs for Cas9, Cpf1 or Custom are provided. For Cas9, the length of protospacer (20nt) and PAM sequence (NGG) are fixed. The length of protospacer and PAM sequence can be set by users for the other two. Simple rules are:
	*  Cas9	:On-target: 20 nt protospacer + NGG, off-target: 20 nt + NRG
	*  Cpf1	:On-target: TTTN/TTN + 23/24/25 nt protospacer
	*  Custom	:On-target: 15-25 nt protospacer + custom PAM sequence
* Screening for sgRNAs that are capable of simultaneously targeting multiple genes.
* Adopting the “one-for-all” strategy to build a whole-genome sgRNA database with high efficiency.
* Running offline, with both command-line and graphical user interface modes, and the ability to export multiple formats for further batch analysis or visualization.
* A substantial pre-stored sgRNA database for 71 public plant reference genomes. CRISPR-Local has been applied to 71 public plant genomes, for both CRISPR/Cas9 (with NGG PAM) and CRISPR/cpf1 (with TTV(A/G/C) and TTTV PAM), which can be directly searched or downloaded from http://crispr.hzau.edu.cn/CRISPR-Local/.

### A General View of CRISPR-Local
* UD-build could build another modified sgRNA database based on the additional sequencing data on non-reference lines, and be called as User’s sgRNA database (UD). Both aligned reads and unmapped reads are used to screen suitable sgRNA, the difference is that each aligned read can be assigned to one gene according to its alignment position. User can summit their own sequence data in bam/sam/fasta/fastq format.
* Program DB-search would compare the RD with the UD, and the results consist of three parts: The sgRNAs present exclusively in RD (RD only, RO) or UD (UD only, UO), or both (BO). Generally, the sgRNAs from BO is preferential, and cautions should be required for UO, and RO is strongly not recommended especially when the supplied NGS data is adequate enough. 
* PL-search is typically designed for editing multiple paralogs to address redundancy concerns, like knocking out multiple key genes for several pathways determining one phenotype. This step follows the function of above DB-search mode that distinguish targets to RO, UO or BO. Moreover, PL-search will additionally partition the targets to common (suits to all candidates) and exclusive (individually matched) according to any submitted list of genes.

### Prerequisites:
The following software and libraries are additionally required: Seqmap (version: 1.0.12), Samtools (>1.3.1), Python (>2.7) with the scikit-learn (0.16.1), biopython, pandas, keras, numpy, and scipy libraries, and Perl (>5.10) with the Parallel::ForkManager, File::Basename, Getopt::Long, Data::Dumper, Cwd modules. Most of the installation steps are fully automatic using a simple command line on a Linux system. 

### 1.Prepare CRISPR-Local input (fasta/genome) files
(1) Reference genome sequences with fasta file, can be downloaded from Ensembl Plants, NCBI or other sources.
  * Before performing CRISPR-Local program, multiple fasta files need be combined into one big fasta file by using cat command.
 
(2) The reference annotation file in GFF3 format, can also be downloaded from Ensembl Plants.
  * The chromosome names (sequence headers) in the FASTA reference genome file must be the same as the first column of GFF3 annotation file.

  * Example GFF3 file
```
##gff-version 3
#!genome-build  Pmarinus_7.0
#!genome-version Pmarinus_7.0
#!genome-date 2011-01
#!genebuild-last-updated 2013-04
GL476399 Pmarinus_7.0 supercontig 1 4695893 . . . ID=supercontig:GL476399;Alias=scaffold_71
GL476399 ensembl gene 2596494 2601138 . + . ID=gene:ENSPMAG00000009070;Name=TRYPA3;biotype=protein_coding;description=Trypsinogen A1%3B Trypsinogen a3%3B Uncharacterized protein [Source:UniProtKB/TrEMBL%3BAcc:O42608];logic_name=ensembl;version=1
GL476399 ensembl transcript 2596494 2601138 . + . ID=transcript:ENSPMAT00000010026;Name=TRYPA3-201;Parent=gene:ENSPMAG00000009070;biotype=protein_coding;version=1
GL476399 ensembl exon 2596494 2596538 . + . Name=ENSPMAE00000087923;Parent=transcript:ENSPMAT00000010026;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;rank=1;version=1
GL476399 ensembl exon 2598202 2598361 . + . Name=ENSPMAE00000087929;Parent=transcript:ENSPMAT00000010026;constitutive=1;ensembl_end_phase=2;ensembl_phase=1;rank=2;version=1
GL476399 ensembl exon 2599023 2599282 . + . Name=ENSPMAE00000087937;Parent=transcript:ENSPMAT00000010026;constitutive=1;ensembl_end_phase=1;ensembl_phase=2;rank=3;version=1
GL476399 ensembl exon 2599814 2599947 . + . Name=ENSPMAE00000087952;Parent=transcript:ENSPMAT00000010026;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;rank=4;version=1
GL476399 ensembl exon 2600895 2601138 . + . Name=ENSPMAE00000087966;Parent=transcript:ENSPMAT00000010026;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;rank=5;version=1
GL476399 ensembl CDS 2596499 2596538 . + 0 ID=CDS:ENSPMAP00000009982;Parent=transcript:ENSPMAT00000010026
GL476399 ensembl CDS 2598202 2598361 . + 2 ID=CDS:ENSPMAP00000009982;Parent=transcript:ENSPMAT00000010026
GL476399 ensembl CDS 2599023 2599282 . + 1 ID=CDS:ENSPMAP00000009982;Parent=transcript:ENSPMAT00000010026
GL476399 ensembl CDS 2599814 2599947 . + 2 ID=CDS:ENSPMAP00000009982;Parent=transcript:ENSPMAT00000010026
GL476399 ensembl CDS 2600895 2601044 . + 0 ID=CDS:ENSPMAP00000009982;Parent=transcript:ENSPMAT00000010026
```
	
### 2. How to run CRISPR-Local

####  (1) program RD-build:
```
Parameters
--- Required ---
	
	-m <string>	:Given specific PAM for sgRNA design : Cas9, Cpf1 or Custom (default: Cas9)

	  Cas9		:On-target: 20 nt protospacer + NGG, off-target: 20 nt + NRG
	  Cpf1		:On-target: TTTN/TTN + 23/24/25 nt protospacer
	  Custom	:On-target: 15-25 nt protospacer + custom PAM sequence

	-i <string>	:reference genome sequence file (In *.fasta format)
	-g <string>	:reference genome annotation file (In *.gff3 format)

--- Options ---

	-o <string>	:Output path (default: current directory)
	-l <string>	:Name prefix for output file (default: CRISPR)
	-p <int>	:The number of process to use (default: 1)

	For Cas9:

	 -U <int>	:An integer specifying the number of base pairs expanding 5'-end of exon boundary (default: 0, >15 is not allowed)
	 -D <int>	:An integer specifying the number of base pairs expanding 3'-end of exon boundary (default: 0, >3 is not allowed)

	For Cpf1 or Custom:

	  -x <int>	:Cpf1: Length of spacer: between 23 to 25 (default: 24 nt);
			:Custom: Length of spacer: between 15 to 25 (default: 20 nt);
	  -t <string>	:Cpf1: Type of PAM sequence: TTX or TTTX (X: One of A C G T R Y M K S W H B V D N; default: TTTN)
			:Custom: Type of PAM sequence (default: NGG)

	-h: show this help
```
##### (i) Designing sgRNA for Cas9:

* In cas9 mode, this script: 
	* (a) Screening all possible on-target sgRNAs (with NGG PAM type) and scoring them with the Rule set 2 algorithm (Doench et al., 2016) for each exon, and screening all potential off-target sites (with NRG PAM type); 
	* (b) Retrieving potential off-target sites against each candidate sgRNA by using SeqMap program (Jiang and Wong, 2008) under a default maximum mismatch of 4; 
	* (c) Appling the cutting frequency determination (CFD) score (Doench et al., 2016) to predict the effects of each off-target site, the highest CFD score for each sgRNA is retained, and all the genome-wide results were exported as the RD.

* In addition, the number of base pairs expanding 5'- and 3'-end of exon boundary can be assigned to indicate if the sgRNA can be accepted when they are partially overlapped with intron while the likely editing positions are still located in exon to knock out likewise.

```
Example:
----------------------------
perl RD-build.pl -m cas9 -i Reference_Genome.fa -g Reference_annotation.gff3 -o /opt/your_dir/ -l Label -U 15 -D 3 -p 8
    
* This command will generate a reference database file and a log file.

Example of reference database files:
----------------------------
	
TAIR10.reference.database.txt

AT1G01010  1:+3855  CGTTGAAGTAGCCATCAGCG+AGG  0.7909  AT1G06430  1:+1962511	CGATGAAGCAGCCATCTGCA+CAG  4  AT1G01010.1.exon1;  [3760:3631:283:224:95];   0.0214
AT1G01010  1:-3849  TGATGGCTACTTCAACGTCG+CGG  0.7724  AT1G01740  1:+274664	CTTTGGCTACTTCAACATCG+CAG  4  AT1G01010.1.exon1;  [3760:3631:283:64:-65];	  0.0932
AT1G01010  1:+4118  GTTGAGGTCAAGGACCAGTG+GGG  0.7487  AT1G51035  1:+18917578	GTTCAGTTCATGAACCAGTG+CAG  4  AT1G01010.1.exon2;  [3760:3996:281:122:358];  0.0223
AT1G01010  1:-5499  TTCACCGTGTTGGTGGATGG+AGG  0.7036  AT1G31080  1:+11092019	GTCACCGTCTTGGTGGATCC+CGG  4  AT1G01010.1.exon6;  [3760:5439:461:400:2079]; 0.0990
AT1G01010  1:+4102  GCTTACCGGAGAATCTGTTG+AGG  0.7031  AT1G04680  1:+1307049	GCTAACCGGAGAAACCGTTA+GAG  4  AT1G01010.1.exon2;  [3760:3996:281:106:342];  0.0478
AT1G01010  1:+3690  CAGAGAGCGAGAGAGATCGA+CGG  0.7005  AT1G30540  1:+10817050	AAGAGAGAGAGAGAGAGAGA+GAG  4  AT1G01010.1.exon1;  [3760:3631:283:59:-70];	  0.0107


* There are 11 columns in the RD file and their corresponding meanings are listed below:

Column 1:	The name of gene where the sgRNA located.
Column 2:	The chromosome and the coordinate of the start position of the sgRNA.
Column 3:	The sequence of sgRNA(23nt).
Column 4:	The on-target score of the sgRNA. 
Column 5:	The name of off-target gene with the highest CFD score.
Column 6:	The chromosome and the coordinate of the start position of the off-target site with the highest CFD score.
Column 7:	The sequence of off-target site.
Column 8:	The number of mismatches between sgRNA and off-target site.
Column 9:	The name of exon where the sgRNA located(split by ;).
Column 10:	The number that split by ":" means "TSS position", "exon start position", "length of exon", "relative position of sgRNA against exon" and "relative position of sgRNA against TSS", respectively.
Column 11:	The highest CFD score between sgRNA and all off-target sites.
```
##### (ii) Cpf1
* Cpf1 is a single RNA-guided endonuclease lacking tracrRNA, and it utilizes a T-rich PAM (5'-TTN or 5'-TTTN). Each mature crRNA begins with 19 nt of the direct repeat followed by 23¨C25 nt of the spacer sequence (Zetsche and Gootenberg et al.,2015);
* In cpf1 mode, this script: 
	* (a) Screening all possible on-target sgRNAs (with TTX or TTTX PAM type, X represents any base including A C G T R Y M K S W H B V D N),and screening all potential off-target sites (with TTX or TTTX PAM type). The type of PAM and length of protospacer are determined by user;
	* (b) Scoring all sgRNAs with the method raised by Kim et al (2017);
	* (c) Retrieving potential off-target sites against each candidate sgRNA by using SeqMap program under a default maximum mismatch of 4; 
	* (d) Summarizing the off-target results, and all the genome-wide results were exported as the RD.
```
Example:
----------------------------
perl RD-build.pl -m cpf1 -i Reference_Genome.fa -g Reference_annotation.gff3 -o /opt/your_dir/ -l Label -x 24 -t TTTV -p 8

* This command will generate a reference database file and a log file.

Example of reference database files:
----------------------------
	
TAIR10.reference.database.txt

AT1G67220	1:-25147036	TTTC+CTCGGTTTGAATCTTTCCTTTGTT	35.156925	5	U1	0,1,4,0,0	NA	AT1G67220.1.exon1;	[25145587:25145587:2181:731:731];	NA
AT1G67220	1:-25147059	TTTC+TTCTCATAGTTCAAAGACCTTTCC	56.483639	5	R1	0,2,3,0,0	NA	AT1G67220.1.exon1;	[25145587:25145587:2181:708:708];	NA
AT1G67220	1:-25147095	TTTC+ATCGGCTCAACAATATCCACACCA	55.160500	8	R0	2,3,2,1,0	AT1G67220(1:-25146972);AT1G67220(1:-25147302);	AT1G67220.1.exon1;	[25145587:25145587:2181:672:672];	NA
AT1G67220	1:-25147164	TTTC+ATTGGCTCAACAATAACCACATCA	48.010044	8	R3	0,0,0,7,1	NA	AT1G67220.1.exon1;	[25145587:25145587:2181:603:603];	NA
AT1G67220	1:-25147173	TTTA+TTACATTTCATTGGCTCAACAATA	42.615669	0	NM	0,0,0,0,0	NA	AT1G67220.1.exon1;	[25145587:25145587:2181:594:594];	NA
AT1G67220	1:-25147233	TTTC+ATTGGCTCAACAATATTCACACCA	41.623928	8	R2	0,0,3,4,1	NA	AT1G67220.1.exon1;	[25145587:25145587:2181:534:534];	NA
AT1G67220	1:-25147302	TTTC+ATCGGCTCAACAATATCCACACCA	55.647942	8	R0	2,3,2,1,0	AT1G67220(1:-25146972);AT1G67220(1:-25147095);	AT1G67220.1.exon1;	[25145587:25145587:2181:465:465];	NA
* There are 11 columns in the RD file and their corresponding meanings are listed below:

Column 1:	The name of gene where the sgRNA located.
Column 2:	The chromosome and the coordinate of the start position of the sgRNA.
Column 3:	The sequence of sgRNA.
Column 4:	The on-target score of the sgRNA. 
Column 5:	The number of off-target sites.
Column 6:	Type of match between sgRNA and off-target sites (NM:no match found; U0:Best match found was a unique exact match; U1:Best match found was a unique 1-error match; U2:Best match found was a unique 2-error match... R0:Multiple exact matches found; R1:Multiple 1-error matches found, no exact matches; R2:Multiple 2-error matches found, no exact or 1-error matches.)
Column 7:	The number of exact, 1-error, 2-error, 3-error and 4-error matches found.
Column 8:	The gene and position in which exact match was found. (If there is no exact match, then denoted by NA)
Column 9:	The name of exon where the sgRNA located(split by ;).
Column 10:	The number that split by ":" means "TSS position", "exon start position", "length of exon", "relative position of sgRNA against exon" and "relative position of sgRNA against TSS", respectively.
Column 11:	The highest off-target score between sgRNA and all off-target sites.(There is no available off-target scoring method for Cpf1 sgRNA, denoted by NA)
```
####  (2) Program UD-build (if necessary):

* This script use to accept user's data to build user's sgRNA database.
* According to the contents and format of the data file, different methods were used to design sgRNAs and build database.
* Program UD-build supports input alignment file in Bam or Sam format, and raw reads file in fasta/fa/fasta.gz/fa.gz or fq/fastq/fq.gz/fastq.gz format.
* If your file is Bam or Sam format, the annotation file in GFF3 format is needed to be specified.

 ```    
Parameters
--- Required ---

	-m <string>	:Given specific PAM for sgRNA design : Cas9, Cpf1 or Custom (default: Cas9)

	  Cas9		:On-target: 20 nt protospacer + NGG, off-target: 20 nt + NRG
	  Cpf1		:On-target: TTTN/TTN + 23/24/25 nt protospacer
	  Custom	:On-target: 15-25 nt protospacer + custom PAM sequence

	-i <string>	:User's sequence file (In bam, fastq or fasta format)
	-g <string>	:Genome annotation file (This parameter is required only if the format of sequence file is bam)

--- Options ---

	-o <string>	:Output path (default: current directory)
	-p <int>	:The number of process to use (default:1)

	-h: show this help
	
	For Cpf1 or Custom mode:

	  -x <int>	:Cpf1: Length of spacer: between 23 to 25 (default: 24 nt);
			:Custom: Length of spacer: between 15 to 25 (default: 20 nt);
	  -t <string>	:Cpf1: Type of PAM sequence: TTX or TTTX (X: One of A C G T R Y M K S W H B V D N; default: TTTN)
			:Custom: Type of PAM sequence (default: NGG)

Example:
----------------------------
For Cas9 mode:

	perl UD-build.pl -m cas9 -i Your_data.bam -g Annotation.gff3 -o /your_dir/ -p 10

For Cpf1 mode:

	perl UD-build.pl -m cpf1 -i Your_data.fasta -o /your_dir/ -p 10 -x 24 -t TTTV

For Custom mode:

	perl UD-build.pl -m custom -i Your_data.fastq -o /your_dir/ -t NRG -p 10 -x 20

* Two user's sgRNA database files would be generated when the input file is Bam or Sam format. 

Example of user's database files:
----------------------------
ZmC01.gene.sgRNA.db.alignment.txt
      
Zm00001d022658  1:+4991138      5       0.2965  ATTCTGATTATATAGATATT+AGG
Zm00001d022658  1:+4991139      1       0.2965  ATTCTGATTATATAGATATT+AGG
Zm00001d022658  1:+4991231      4       0.5081  TGTTAGCCATGACATGTTTG+AGG
Zm00001d022658  1:+4991232      4       0.5198  GTTAGCCATGACATGTTTGA+GGG
Zm00001d022658  1:+4991233      4       0.6922  TTAGCCATGACATGTTTGAG+GGG

ZmC01.intergenic.sgRNA.db.alignment.txt

Intergenic      1:+1024980      1       0.2012  CAGTCGTTGCCAAGCGTTCT+TGG
Intergenic      1:+1025011      2       0.4447  ACCAGCAAGCAGCGCACCAC+CGG
Intergenic      1:+1025033      2       0.6557  GCAAGCAGCGCACCACCACA+AGG
Intergenic      1:+1025039      2       0.5222  AGCGCACCACCACAAGGTTG+CGG
Intergenic      1:+1025057      2       0.5021  TGCGGCTTCGAGCACTGCAC+CGG
Intergenic      1:+1025080      1       0.5388  CAAGCAGCGCACCACCAGCA+AGG
Intergenic      1:+1025085      1       0.5927  AGCGCACCACCAGCAAGGAG+TGG


* There are 5 columns in the UD file and their meanings are listed below:

Column 1:	The name of gene where the sgRNA located ("Intergenic" means sgRNA locate in intergenic region).
Column 2:	The chromosome and the coordinate of the start position of the sgRNA.
Column 3:	The number of reads that contains this sgRNA.
Column 4:	The on-target score of the sgRNA.
Column 5:	The sequence of sgRNA (23nt).

* A user's sgRNA database file would be generated when the input file is fasta or fastq format. 

Example of user's database files:
----------------------------
ZmC01.gene.sgRNA.db.fastq.txt 

ZmC01        3       0.5809  ACAAACAGAGGTCTAAAGCA+AGG
ZmC01        5       0.4756  ACAAACAGCCGGTGAAGCTC+CGG
ZmC01        1       0.4811  ACAAACATTACCTTGTTGAG+AGG
ZmC01        1       0.5456  ACAAACCTGCTCTCAGGGGT+GGG
ZmC01        1       0.5662  ACAAACCTTTCTGTTCTGAT+GGG
ZmC01        2       0.7371  ACAAACGCATGATACATAGG+TGG
ZmC01        16      0.4073  ACAAACGGCCGGCGGCAGCT+AGG

Column 1:	The ID of sgRNA.
Column 2:	The number of reads that contains this sgRNA.
Column 3:	The on-target score of the sgRNA. 
Column 4:	The sequence of sgRNA(23nt).

ZmC01.gene.sgRNA.db.fasta.txt
        
chr1    +1124   0.1066  TAATCAAATAAATAAGTTTA+TGG
chr1    +1279   0.3327  AGTAATACATTCTTATAAAA+TGG
chr1    +1428   0.4202  GAGTCAGTGTCGTTATGTTA+TGG
chr1    +1501   0.4973  TTACAAGGGAAGTCCCCAAT+TGG
chr1    +1568   0.2559  AATCTTCTAATTACTGTATA+TGG
chr1    +1620   0.3354  GTGGCCAAGGTTCCGTCATT+TGG
chr1    +1699   0.3337  ACATCTATCTCCATATGATA+TGG

Column 1:	The ID of reads.
Column 2:	The coordinate of the start position of the sgRNA.
Column 3:	The on-target score of the sgRNA. 
Column 4:	The sequence of sgRNA(23nt).
```
#### (3) Program DB-search:

* User could provide a gene list file in TXT format, one gene per line, the program DB-search would search and compare RD and UD.
* If you did not run the program UD-build, you can specify only the RD.
* DB-search select the important columns in the database file and finally output the result.
```
Parameters
--- Required ---

	-g <string>	:Query gene list file
	-i <string>	:Reference sgRNA database (RD)

--- Optional ---

	-u <string>	:User's sgRNA database (UD)
	-o <string>	:Output path (defualt: current directory)
	-l <label>	:Name prefix for output file (default:DB_search)

	-h		:show this help

Example:
     ----------------------------
     perl DB-search.pl -g ZmB73_query_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/ -l label

     Output of this command
     Invalid_gene_RD.list	(Includes genes that do not exist in RD)
     label_result_RO.txt	(The RD-specific sgRNAs)
     label_result_UO.txt	(The UD-specific sgRNAs)
     label_result_BO.txt	(The common sgRNAs)

     Zm00001d001775  2:-547319       0.6568  0.0728  CGGGCGCATCATGCGCCGCG+CGG
     Zm00001d001775  2:-547292       0.6551  0.0458  GGAGAACGGAAAGATCGCTA+GGG
     Zm00001d001775  2:+547304       0.6542  0.0543  TTCCGTTCTCCGTCACCGCG+CGG
     Zm00001d001792  2:+1029153      0.7932  0.0457  GCCGGAGTACTCGAGCAGCG+CGG
     Zm00001d001792  2:+1027126      0.7213  0.1346  GAAAGTGAGATACAAGCCAG+TGG
     Zm00001d001792  2:+1028960      0.7101  0.0478  CGGCAGGTGATGAGTCCTCG+GGG
     
    * User could select the sgRNA with high on-target score (column 3) and low off-target score (column 4).
```   
#### (4) Program PL-search:
    
* PL-search is a local tool for search exclusive and common target for paralogous gene pair.
* Example of paralog gene list:
```
Zm00001d049540,Zm00001d024543
Zm00001d048890,Zm00001d029176
Zm00001d031930,Zm00001d018487,Zm00001d003312
Zm00001d036877,Zm00001d046535
```
```     
Parameters
--- Required ---

	-g <string>	:Paralog gene list
	-i <string>	:Reference sgRNA database (RD)

--- Optional ---

	-u <string>	:User's sgRNA database (UD)
	-o <string>	:Output path (defualt: current directory)
	-l <label>	:Name prefix for output file (default:PL_search)

	-h: show this help and exit

Example:
     ----------------------------
      perl PL-search -g ZmB73_paralog_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/ -l label
      
      Output of this command
      label_common_targets.txt
      label_exclusive_targets.txt
      
      label_common_targets.txt
      
      GAAAATGTTGCCCATCGATA+TGG UD      Zm00001d031930:1:-207131656:0.428795    Zm00001d018487:5:-221689986:0.428795    Zm00001d003312:2:+39761372:0.428795
      AAGAAAGGGCTGCCCATTCT+TGG UD      Zm00001d031930:1:-207131394:0.351515    Zm00001d018487:5:-221689268:0.351515    Zm00001d003312:2:+39761819:0.351515
      CTCCTCCGGCAGCGGGAGCT+GGG UD      Zm00001d036877:6:+105440097:0.354041    Zm00001d046535:9:-94504945:0.354041
      CGCCACCCAGCTCCCGCTGC+CGG UD      Zm00001d036877:6:-105440102:0.340637    Zm00001d046535:9:+94504940:0.340637
      TGCAGGCGGCGTACATCCTG+TGG UD      Zm00001d036877:6:-105440202:0.487230    Zm00001d046535:9:+94504840:0.487230
      CCCGCTGCCGGAGGAGTTCC+TGG UD      Zm00001d036877:6:-105440090:0.158463    Zm00001d046535:9:+94504952:0.163690
      
      Column 1: The sequence of sgRNA.
      Column 2:	The sgRNAs exist in UD are indicated with "UD" mark.
      Column 3,4,5:	The gene ID; the coordinate o sgRNA; the on-target score of the sgRNA. 
      
      label_exclusive_targets.txt
      
      Zm00001d031930  1:-207131270    0.230513        TGGAGTGGATCCAGCTATTA+TGG Zm00001d028694  1:+43348332     TGTAGTTGATCCAGGTACTA+CAG 4       0.0016  RD      0
      Zm00001d031930  1:+207131844    0.139577        TTCCCAAGCAATGGAATTTA+TGG Zm00001d024478  10:-73139286    TTTACAAGCAGTGGAACTTA+AGG 4       0.2656  UD      1
      Zm00001d018487  5:+221690821    0.782208        GAACCTTGAGAACAACAGTG+AGG Zm00001d027293  1:+2033102      TTCCCTTGACAACAACAGTG+CGG 4       0.1247  RD      0
      Zm00001d018487  5:-221690848    0.776871        CTTCAAGGACACCTACCCAG+AGG Zm00001d023582  10:+10766178    TTTCAAGCACCCCTTCCCAG+TGG 4       0.0492  UD      312
      Zm00001d018487  5:+221688910    0.748847        GCAGCAATAGACAAACTGAG+AGG Zm00001d024647  10:+82039457    GAAGCTATAGACAAGCTCAG+GAG 4       0.0417  UD      303
      Zm00001d018487  5:+221690357    0.733564        GTAGCCCAGATGAACACTGG+CGG Zm00001d025888  10:+133099566   GTAGGCAAAATGAACTCTGG+TAG 4       0.0000  UD      357
      
      Column 1-11: As stated above.
      Column 12 : The sgRNAs exist in UD are indicated with "UD" mark, else indicated with "RD" mark.
      Column 13 : The frequency of sgRNA in UD.
```

For other technical issues about CRISPR-P 2.0, please contact: liuhao1122@webmail.hzau.edu.cn
