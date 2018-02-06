# CRISPR-Local

## A local tool for high-throughput CRISPR single-guide RNA (sgRNA) design in plants.

The following additional software and libraries are required: Seqmap (version: 1.0.12), Samtools (>1.3.1), Python (>2.7) with the scikit-learn (0.16.1), biopython, pandas, numpy, and scipy libraries, and Perl (>5.10) with the Parallel::ForkManager, File::Basename, Getopt::Long, Data::Dumper, Cwd modules. Most of the installation steps are fully automatic using a simple command line on a Linux system. 

### Highlight:

* In the plants for gene editing, the genome of the material and the reference genome are quite different. Leading to the fact that the sgRNA designed according to the reference genome is not well used in other materials.
* User can summit their own sequence data in bam/sam/fasta/fastq format, then use program UD-build to make user's sgRNA database (UD).
* Program DB-search would compare the RD with the UD, and output the RD-, UD-specific sgRNAs and commom sgRNAs, respectively. In general, the result of UD-specific sgRNAs and Common sgRNAs is better than RD-specific sgRNA.
* PL-search mode keeps the function of DB-search mode that design targets in RD or both RD and UD. Moreover, PL-search would design targets for paralogs. In this mode, paralogous gene list is required, then users can design common targets or exclusive targets according to their needs.

### 1.Prepare CRISPR-Local input (fasta/genome) files

(1) Reference genome fasta file, can be downloaded from Ensembl Plants or NCBI or other source.
  * Before perform CRISPR-Local program, multiple fasta files need be combined into one big fasta file by using cat command.
  * In order for this step to work correctly, the chromosome names (sequence headers) in the FASTA reference genome file must be the same as the first column of GFF3 annotation file.
(2) The reference annotation file in GFF3 format can be downloaded from Ensembl Plants.
  
  * Example GFF3 file
```r	
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

####  (1) progam RD-build:

* This script use to find all potential off-target sites from every gene sequence, with NGG and NAG PAM type(20mer(NGG/NAG)), and all possible sgRNAs with NGG PAM from every exon sequence(NNNN20merNGGNNN).
* RD-build score all the sgRNAs by Rule set2 algorithm(John G Doench et al. 2016), then call the SeqMap(Jiang et al. 2008) program to identify how much potential off-target site does each sgRNA have, with maximum number of mismatches up to 4. 
* CFD score(John G Doench et al. 2016) was used to predict sgRNA off-target effects of each sgRNA and its potential off-target site, keeping one result with highest CFD score of each sgRNA, then outputting the results to the reference sgRNA database (RD)
* In addition, user can specify the number of bases to expanding 5'-end and 3'-end for each exon respectively.
```r
Example:
----------------------------
perl RD-build.pl -i Reference_Genome.fa -g Reference_annotation.gff3 -o /opt/your_dir/ -l Label -U 15 -D 3 -p 8
    
* This command will generate a reference database file and a log file.

Example of reference database files:
----------------------------
	
TAIR10.reference.database.txt

AT1G01010  1:+3855  CGTTGAAGTAGCCATCAGCGAGG  0.7909  AT1G06430  1:+1962511	CGATGAAGCAGCCATCTGCACAG  4  AT1G01010.1.exon1;  [3760:3631:283:224:95];   0.0214
AT1G01010  1:-3849  TGATGGCTACTTCAACGTCGCGG  0.7724  AT1G01740  1:+274664	CTTTGGCTACTTCAACATCGCAG  4  AT1G01010.1.exon1;  [3760:3631:283:64:-65];	  0.0932
AT1G01010  1:+4118  GTTGAGGTCAAGGACCAGTGGGG  0.7487  AT1G51035  1:+18917578	GTTCAGTTCATGAACCAGTGCAG  4  AT1G01010.1.exon2;  [3760:3996:281:122:358];  0.0223
AT1G01010  1:-5499  TTCACCGTGTTGGTGGATGGAGG  0.7036  AT1G31080  1:+11092019	GTCACCGTCTTGGTGGATCCCGG  4  AT1G01010.1.exon6;  [3760:5439:461:400:2079]; 0.0990
AT1G01010  1:+4102  GCTTACCGGAGAATCTGTTGAGG  0.7031  AT1G04680  1:+1307049	GCTAACCGGAGAAACCGTTAGAG  4  AT1G01010.1.exon2;  [3760:3996:281:106:342];  0.0478
AT1G01010  1:+3690  CAGAGAGCGAGAGAGATCGACGG  0.7005  AT1G30540  1:+10817050	AAGAGAGAGAGAGAGAGAGAGAG  4  AT1G01010.1.exon1;  [3760:3631:283:59:-70];	  0.0107


* There are 11 columns in the RD file and their meanings are listed below:

Column 1:	The name of gene where the sgRNA located.
Column 2:	The chromosome and the coordinate of the start position of the sgRNA.
Column 3:	The sequence of sgRNA(23nt).
Column 4:	The on-target score of the sgRNA. 
Column 5:	The name of off-target gene with the highest CFD score.
Column 6:	The chromosome and the coordinate of the start position of the off-target site with the highest CFD score.
Column 7:	The sequence of off-target site.
Column 8:	The number of mismatches between sgRNA and off-target site.
Column 9:	The name of exon where the sgRNA located(split by ;).
Column 10:	The number that split by ":" means "TSS position", "exon start position", "length of exon", "relative positon of sgRNA against exon" and "relative positon of sgRNA against TSS", respectively.
Column 11:	The highest CFD score between sgRNA and all off-target sites.
```
####  (2) Program UD-build (if nessesary):

* This script use to accept user's data to build user's sgRNA database.
* According to the contents and format of the data file, different methods was used to design sgRNAs and build database.
* Program UD-build supports input alignment file in Bam or Sam format, and raw reads file in fasta/fa/fasta.gz/fa.gz or fq/fastq/fq.gz/fastq.gz format.
* If your file is Bam or Sam format, the annotation file in GFF3 format is needed to be specified.
* Program RD-build and UD-build could be run at the same time.
 ```    
	  Example:
          ----------------------------
	  perl UD-build.pl -i Your_data.bam -g Annotation.gff3 -o /your_dir/ -p 10

	  or 

	  perl UD-build.pl -i Your_data.fasta.gz -o /your_dir/ -p 10

    * Two user's sgRNA database files would be generated when the input file is Bam or Sam format. 

          Example of user's database files:
          ----------------------------
          ZmC01.gene.sgRNA.db.alignment.txt
        
          Zm00001d022658  1:+4991138      5       0.2965  ATTCTGATTATATAGATATTAGG
          Zm00001d022658  1:+4991139      1       0.2965  ATTCTGATTATATAGATATTAGG
          Zm00001d022658  1:+4991231      4       0.5081  TGTTAGCCATGACATGTTTGAGG
          Zm00001d022658  1:+4991232      4       0.5198  GTTAGCCATGACATGTTTGAGGG
          Zm00001d022658  1:+4991233      4       0.6922  TTAGCCATGACATGTTTGAGGGG

          ZmC01.intergenic.sgRNA.db.alignment.txt

	Intergenic      1:+1024980      1       0.2012  CAGTCGTTGCCAAGCGTTCTTGG
        Intergenic      1:+1025011      2       0.4447  ACCAGCAAGCAGCGCACCACCGG
        Intergenic      1:+1025033      2       0.6557  GCAAGCAGCGCACCACCACAAGG
        Intergenic      1:+1025039      2       0.5222  AGCGCACCACCACAAGGTTGCGG
        Intergenic      1:+1025057      2       0.5021  TGCGGCTTCGAGCACTGCACCGG
        Intergenic      1:+1025080      1       0.5388  CAAGCAGCGCACCACCAGCAAGG
        Intergenic      1:+1025085      1       0.5927  AGCGCACCACCAGCAAGGAGTGG

    * There are 5 columns in the UD file and their meanings are listed below:

        Column 1:	The name of gene where the sgRNA located("Intergenic" means sgRNA locate in intergenic region).
        Column 2:	The chromosome and the coordinate of the start position of the sgRNA.
        Column 3:	The number of reads that contains this sgRNA.
        Column 4:	The on-target score of the sgRNA. 
        Column 5:	The sequence of sgRNA(23nt).

    * A user's sgRNA database file would be generated when the input file is fasta or fastq format. 

        Example of user's database files:
        ----------------------------
        ZmC01.gene.sgRNA.db.fastq.txt 

        ZmC01_sgRNA        3       0.5809  ACAAACAGAGGTCTAAAGCAAGG
        ZmC01_sgRNA        5       0.4756  ACAAACAGCCGGTGAAGCTCCGG
        ZmC01_sgRNA        1       0.4811  ACAAACATTACCTTGTTGAGAGG
        ZmC01_sgRNA        1       0.5456  ACAAACCTGCTCTCAGGGGTGGG
        ZmC01_sgRNA        1       0.5662  ACAAACCTTTCTGTTCTGATGGG
        ZmC01_sgRNA        2       0.7371  ACAAACGCATGATACATAGGTGG
        ZmC01_sgRNA        16      0.4073  ACAAACGGCCGGCGGCAGCTAGG

        Column 1:	The ID of sgRNA.
        Column 2:	The number of reads that contains this sgRNA.
        Column 3:	The on-target score of the sgRNA. 
        Column 4:	The sequence of sgRNA(23nt).

	ZmC01.gene.sgRNA.db.fasta.txt
        
	chr1    +1124   0.1066  TAATCAAATAAATAAGTTTATGG
        chr1    +1279   0.3327  AGTAATACATTCTTATAAAATGG
        chr1    +1368   0.3355  TAAATGAGATGTTGAATTAGAGG
        chr1    +1428   0.4202  GAGTCAGTGTCGTTATGTTATGG
        chr1    +1501   0.4973  TTACAAGGGAAGTCCCCAATTGG
        chr1    +1568   0.2559  AATCTTCTAATTACTGTATATGG
        chr1    +1620   0.3354  GTGGCCAAGGTTCCGTCATTTGG
        chr1    +1699   0.3337  ACATCTATCTCCATATGATATGG

        Column 1:	The ID of reads.
        Column 2:	The coordinate of the start position of the sgRNA.
        Column 3:	The on-target score of the sgRNA. 
        Column 4:	The sequence of sgRNA(23nt).
```
#### (3) Program DB-search:

* User could provide a gene list file in TXT format, one gene per line, the program DB-search would search and compare RD and UD.
* If you did not run the program UD-build, you can specify only the RD.
* DB-search select the important columns in the database file, sorted by score, and finally output the result.
```
     Example:
     ----------------------------
     perl DB-search.pl -l ZmB73_query_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/ -N 3

     Output of this command
     Invalid_gene_RD.list	(Includes genes that do not exist in RD)
     Gene_search_result_RD_only.txt	(The RD-specific sgRNAs)
     Gene_search_result_UD_only.txt	(The UD-specific sgRNAs)
     Gene_search_result_Co-sgRNA.txt	(The commom sgRNAs)

     Zm00001d001775  2:-547319       0.6568  0.0728  CGGGCGCATCATGCGCCGCGCGG
     Zm00001d001775  2:-547292       0.6551  0.0458  GGAGAACGGAAAGATCGCTAGGG
     Zm00001d001775  2:+547304       0.6542  0.0543  TTCCGTTCTCCGTCACCGCGCGG
     Zm00001d001792  2:+1029153      0.7932  0.0457  GCCGGAGTACTCGAGCAGCGCGG
     Zm00001d001792  2:+1027126      0.7213  0.1346  GAAAGTGAGATACAAGCCAGTGG
     Zm00001d001792  2:+1028960      0.7101  0.0478  CGGCAGGTGATGAGTCCTCGGGG
     
    * User could select the sgRNA with high on-target score (column 3) and low off-target score (column 4).
```   
#### (4) Program PL-search:
    
* PL-search is a local tool for search exclusive and commom target for paralogous gene pair.
```     
      Example:
     ----------------------------
      perl PL-search -l ZmB73_paralogous_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/
      
      Output of this command
      Paralog_search_result_common.txt
      Paralog_search_result_exclusive.txt
      
      Paralog_search_result_common.txt
      
      GAAAATGTTGCCCATCGATATGG UD      Zm00001d031930:1:-207131656:0.428795    Zm00001d018487:5:-221689986:0.428795    Zm00001d003312:2:+39761372:0.428795
      AAGAAAGGGCTGCCCATTCTTGG UD      Zm00001d031930:1:-207131394:0.351515    Zm00001d018487:5:-221689268:0.351515    Zm00001d003312:2:+39761819:0.351515
      CTCCTCCGGCAGCGGGAGCTGGG UD      Zm00001d036877:6:+105440097:0.354041    Zm00001d046535:9:-94504945:0.354041
      CGCCACCCAGCTCCCGCTGCCGG UD      Zm00001d036877:6:-105440102:0.340637    Zm00001d046535:9:+94504940:0.340637
      TGCAGGCGGCGTACATCCTGTGG UD      Zm00001d036877:6:-105440202:0.487230    Zm00001d046535:9:+94504840:0.487230
      CCCGCTGCCGGAGGAGTTCCTGG UD      Zm00001d036877:6:-105440090:0.158463    Zm00001d046535:9:+94504952:0.163690
      
      Column 1: The sequence of sgRNA(23nt).
      Column 2:	The sgRNAs exist in UD are indicated with "UD" mark.
      Column 3,4,5:	The gene ID; the coordinate o sgRNA; the on-target score of the sgRNA. 
      
      Paralog_search_result_exclusive.txt
      
      Zm00001d031930  1:-207131270    0.230513        TGGAGTGGATCCAGCTATTATGG Zm00001d028694  1:+43348332     TGTAGTTGATCCAGGTACTACAG 4       0.0016  RD      0
      Zm00001d031930  1:+207131844    0.139577        TTCCCAAGCAATGGAATTTATGG Zm00001d024478  10:-73139286    TTTACAAGCAGTGGAACTTAAGG 4       0.2656  UD      1
      Zm00001d018487  5:+221690821    0.782208        GAACCTTGAGAACAACAGTGAGG Zm00001d027293  1:+2033102      TTCCCTTGACAACAACAGTGCGG 4       0.1247  RD      0
      Zm00001d018487  5:-221690848    0.776871        CTTCAAGGACACCTACCCAGAGG Zm00001d023582  10:+10766178    TTTCAAGCACCCCTTCCCAGTGG 4       0.0492  UD      312
      Zm00001d018487  5:+221688910    0.748847        GCAGCAATAGACAAACTGAGAGG Zm00001d024647  10:+82039457    GAAGCTATAGACAAGCTCAGGAG 4       0.0417  UD      303
      Zm00001d018487  5:+221690357    0.733564        GTAGCCCAGATGAACACTGGCGG Zm00001d025888  10:+133099566   GTAGGCAAAATGAACTCTGGTAG 4       0.0000  UD      357
      
      Column 1-11: As stated above.
      Column 12 : The sgRNAs exist in UD are indicated with "UD" mark, else indicated with "RD" mark.
      Column 13 : The frequency of sgRNA in UD.
```
Please forward any question and suggestion about CRISPR-Local to: sunjm0824@webmail.hzau.edu.cn
