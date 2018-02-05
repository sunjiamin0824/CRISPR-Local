#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;
use Data::Dumper;
use Cwd;

local $SIG{__WARN__} = sub {
	my $message = shift;
	die $message;
};

my $time = time();

my($Input,$outdir,$gff3,$process,$opt_help);
GetOptions(
	"i=s"=>\$Input,
	"o=s"=>\$outdir,
	"g:s"=>\$gff3,
	"p:i"=>\$process,
	"h!" =>\$opt_help,
);
my $dir_default = getcwd;             #default output
$outdir ||= $dir_default;
$process ||= "1";

my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,
	perl $0 -i Your_data.bam -g Annotation.gff3 -o /your_dir/ -p 10

	or

	perl $0 -i Your_data.fasta -o /your_dir/ -p 10

--- Options ---

	-h: show this help

	-i:	<User_sequencing_data> (split by ',')
	-g:	<reference_genome_annotation_file>
	-o:	<output_path> (defualt: $dir_default)
	-p:	<string> The number of process to use (default:1)
};

err( $helpmsg ) if $opt_help;

print "\n  Welcome to CRISPR-Local\n";
print "  ---a local tool for high-throughput CRISPR single-guide RNA (sgRNA) design in plants.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 1.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjiamin0824\@qq.com"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today : $local_time\n\n";
print "# Program UD-bulid: CRISPR sgRNA design by using user's data.\n\n";

if (!$Input) {
	print
		"\n" .
		"Please input the option, and press <Enter>.\n" .
		"For quick help, input '-h' and press <Enter>.\n" .
		"\n";
	die "\n";
}

my $pm = Parallel::ForkManager->new($process);

my @input = split /,/,$Input;
foreach my $input (@input) {
	my ($filename,$prefix) = (split /\./,basename($input))[0,-1];
	open  (LOG, ">$outdir/$filename.Log.txt") || die "Can't open $filename.Log.txt for writing!" ."\n";

	print  LOG "######################################### Log #########################################". "\n\n";
	print  LOG "#                                     CRISPR-Local                                        " ."\n";
	print  LOG "#  ---a local tool for high-throughput CRISPR single-guide RNA (sgRNA) design in plants." ."\n";          
	print  LOG "#                                                                                         " ."\n";
	print  LOG "#             contact:  Sun Jiamin, MaizeGo, Email: sunjiamin0824\@qq.com                    \n\n";
	print  LOG "# Time, begin at $local_time."."\n";
	print  LOG "# Program UD-build: CRISPR sgRNA design by using user's data.\n\n";

	system("mkdir -p -m 755 $outdir/$filename.User.sgRNA");

	#############################################################################################################
	print LOG  "# Begin to parse your data...\n";
	print  "# Begin to parse your data...\n";

	if ($prefix =~ /sam|bam/i) {
		print LOG  "# Your data $input is $prefix format.\n";
		print  "# Your data $input is $prefix format.\n";
		if ($prefix =~ /sam/) {
			system "samtools sort -@ $process -o $outdir/$filename.bam $input";
			system "samtools index $outdir/$filename.bam";
			$input = "$outdir/$filename.bam";
		}
		if ($prefix =~ /bam/) {
			unless(-e "$input.bai"){system "samtools index $input";}
		}
		open HEAD, "samtools view -H $input |" or die $!; my @SQ;
		while(<HEAD>) {  
			if(/^\@SQ/) {  
				my ($sq) = $_ =~ /SN:(\S+)/;  
				push (@SQ,$sq);  
				next;
			}
		}
		close HEAD;
		foreach my $chr (@SQ) {
			unless(-e "$outdir/$filename.$chr.sam"){system "samtools view -F 4 $input $chr > $outdir/$filename.$chr.sam";}
			my $line = `wc -l $outdir/$filename.$chr.sam`;
			if ($line =~ /^(\d+)/){
				$line = $1;
			}
			my $file_line = int($line/$process) + 1;
			system("split -d -l $file_line $outdir/$filename.$chr.sam $outdir/$filename.$chr.sam.part");
			unlink("$outdir/$filename.$chr.sam") or die "Can't delete $filename.$chr.sam file";
		}
		##############################################################################################################
		print  "# Reading reference gff3 file.\n";
		print LOG  "# Reading reference gff3 file.\n";
		open (REF,$gff3) or die;
		my (%exact,%left,%right);
		while (<REF>) {
			next if ($_=~/^\#/);
			chomp $_;
			my($chr,$type,$S1,$E1,$info) = (split("\t",$_))[0,2,3,4,8];
			if ($type =~ /gene$/) {
				$chr =~ s/chr//;
				my($attrs) = split( ";", $info );
				$attrs =~ s/id\=//i;
				$attrs =~ s/gene://i;
				my $up = int($S1/1000);
				my $down = int($E1/1000);
				my $mod_up = $S1 % 1000;
				my $mod_down = $E1 % 1000;
				$left{"$chr:$up"} = "$attrs\t$mod_up";
				$right{"$chr:$down"} = "$attrs\t$mod_down";
				next unless ($down - $up)>1;
				foreach ($up+1..$down-1) {
					$exact{"$chr:$_"} = $attrs;
				}
			}
		}
		close REF;
#############################################################################################################
		print LOG  "# Identifying which gene does the reads belongs to.\n";
		print  "# Identifying which gene does the reads belongs to.\n";
		foreach my $chr (@SQ) {
			for (my $P = 0;$P < $process;$P++) {
				$pm ->start and next;
				$P = sprintf("%02d",$P);
				if (-e "$outdir/$filename.$chr.sam.part$P") {
					open(IN,"$outdir/$filename.$chr.sam.part$P") or die;
					open (OUT,">$outdir/$filename.$chr.gene.txt.part$P") or die "Can't open $filename.$chr.gene.txt.part$P for writing!\n";
					open (OUT1,">$outdir/$filename.$chr.intergenic.txt.part$P") or die "Can't open $filename.$chr.intergenic.txt.part$P for writing!\n";
					while (<IN>) {
						chomp;
						my ($str,$chr,$pos,$cigar,$seq)=(split /\t/,$_)[1,2,3,5,9];
						$chr =~ s/chr//;
						if ($cigar =~ /^(\d+)S/) {
							my $clip = $1;
							$pos = $pos - $clip;
						}
						my $n = int($pos/1000);
						my $mod_n = $pos % 1000;
						if (exists $exact{"$chr:$n"}) {
							print OUT $exact{"$chr:$n"}."\t$chr\t$pos\t$seq\n";
						}elsif (exists $left{"$chr:$n"}) {
							my ($chr_up,$up) = (split /\t/,$left{"$chr:$n"})[0,1];
							if ($up < $mod_n) {
								print OUT "$chr_up\t$chr\t$pos\t$seq\n";
							}
						}elsif (exists $right{"$chr:$n"}) {
							my ($chr_down,$down) = (split /\t/,$right{"$chr:$n"})[0,1];
							if ($down > $mod_n) {
								print OUT "$chr_down\t$chr\t$pos\t$seq\n";
							}
						}else{
							print OUT1 "Intergenic\t$chr\t$pos\t$seq\n";
						}
					}
					close IN;
					close OUT;
					close OUT1;
					unlink("$outdir/$filename.$chr.sam.part$P") or die "Can't delete $filename.$chr.sam.part$P file";
				}
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		print LOG  "# Extracting possible sgRNA.\n";
		print  "# Extracting possible sgRNA.\n";
		my $A_pos;my $Reads;my $PAM;my %hash;
		foreach my $chr (@SQ) {
			for (my $P = 0;$P < $process;$P++) {
				$pm->start and next;
				$P = sprintf("%02d",$P);
				if (-e "$outdir/$filename.$chr.gene.txt.part$P") {
					open (IN,"$outdir/$filename.$chr.gene.txt.part$P") or die;
					open (OUT,">$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.fasta.part$P") or die "Can't open $filename.$chr.gene.sgRNA.fasta for writing!\n";
					while (<IN>) {
						chomp;
						my ($gene,$chr,$pos,$seq) = split /\t/,$_;
						while ($seq =~ /(?=(\w{25}GG\w{3}))/g) {
							$A_pos = int($pos + (pos($seq) + 4));
							$Reads = $1;
							$PAM = substr($1,24,3);
							unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
								$hash{"$gene\#$chr:+$A_pos\#$PAM"}{$Reads}++;
							}
						}
						$seq = reverse $seq;
						$seq =~ tr/AGCT/TCGA/;
						while ($seq =~ /(?=(\w{25}GG\w{3}))/g) {
							$A_pos = int($pos + (pos($seq) + 4));
							$Reads = $1;
							$PAM = substr($1,24,3);
							unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
								$hash{"$gene\#$chr:-$A_pos\#$PAM"}{$Reads}++;
							}
						}
					}	
					close IN;
					foreach my $key1 (sort keys %hash) {
						foreach my $key2 (sort keys %{$hash{$key1}}) {
							print OUT ">$key1\#$hash{$key1}{$key2}\n$key2\n";
						}
					}
					close OUT;undef %hash;
					unlink("$outdir/$filename.$chr.gene.txt.part$P") or die "Can't delete $filename.$chr.gene.txt.part$P file";
				}
				if (-e "$outdir/$filename.$chr.intergenic.txt.part$P") {
					open (IN,"$outdir/$filename.$chr.intergenic.txt.part$P") or die;
					open (OUT,">$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.fasta.part$P") or die "Can't open $filename.$chr.intergenic.sgRNA.fasta for writing!\n";
					while (<IN>) {
						chomp;
						my ($chr,$pos,$seq) = (split /\t/,$_)[1,2,3];
						while ($seq =~ /(?=(\w{25}GG\w{3}))/g) {
							$A_pos = int($pos + (pos($seq) + 4));
							$Reads = $1;
							$PAM = substr($1,24,3);
							unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
								$hash{"Intergenic#$chr:+$A_pos\#$PAM"}{$Reads}++;
							}
						}
						$seq = reverse $seq;
						$seq =~ tr/AGCT/TCGA/;
						while ($seq =~ /(?=(\w{25}GG\w{3}))/g) {
							$A_pos = int($pos + (pos($seq) + 4));
							$Reads = $1;
							$PAM = substr($1,24,3);
							unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
								$hash{"Intergenic#$chr:-$A_pos\#$PAM"}{$Reads}++;
							}
						}
					}
					close IN;
					foreach my $key1 (sort keys %hash) {
						foreach my $key2 (sort keys %{$hash{$key1}}) {
							print OUT ">$key1\#$hash{$key1}{$key2}\n$key2\n";
						}
					}
					close OUT;undef %hash;
					unlink("$outdir/$filename.$chr.intergenic.txt.part$P") or die "Can't delete $filename.$chr.intergenic.txt.part$P file";
				}
				$pm->finish;
			}
			$pm->wait_all_children;
		}
				
		my %gene_file;my %inter_file;
		foreach my $chr (@SQ) {
			for (my $P=0;$P < $process;$P++) {
				$P = sprintf("%02d",$P);
				if (-e "$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.fasta.part$P") {
					$gene_file{"$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.fasta.part$P"} = "$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.score.fasta.part$P";
				}
				if (-e "$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.fasta.part$P") {
					$inter_file{"$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.fasta.part$P"} = "$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.score.fasta.part$P";
				}			
			}
		}
		print LOG  "# Calculates the Rule set 2 score for the sgRNA.\n";
		print  "# Calculates the Rule set 2 score for the sgRNA.\n";

		foreach my $key1 (sort keys %gene_file) {
			$pm->start and next;
			my @gene_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$key1","--output","$gene_file{$key1}");
			system(@gene_rs2);
			if($? == -1) {
				die "system @gene_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		foreach my $key2 (sort keys %inter_file) {
			$pm->start and next;
			my @intergene_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$key2","--output","$inter_file{$key2}");
			system(@intergene_rs2);
			if($? == -1) {
				die "system @intergene_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		
		print LOG  "# Output the user's database.\n";
		print  "# Output the user's database.\n";
		open (REG,">$outdir/$filename.User.sgRNA/$filename.gene.sgRNA.db.alignment.txt") or die;
		open (REI,">$outdir/$filename.User.sgRNA/$filename.intergenic.sgRNA.db.alignment.txt") or die;
		foreach my $chr (@SQ) {
			for (my $P = 0;$P < $process;$P++) {
				$P = sprintf("%02d",$P);
				my $pam;
				if (-e "$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.score.fasta.part$P") {
					open (IN, "$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.score.fasta.part$P") or die;
					while (<IN>) {
						chomp;
						if ($_=~s/^>//) {
							my ($id,$pos,$num_score) = (split "#",$_)[0,1,3];
							$pam = (split "#",$_)[2];
							my ($num,$score) = split "_",$num_score;
							print REG "$id\t$pos\t$num\t$score\t";
						}else{
							print REG "$_"."$pam\n";
						}
					}
					close IN;
					unlink ("$outdir/$filename.User.sgRNA/$filename.$chr.gene.sgRNA.score.fasta.part$P") or die;
				}
				if (-e "$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.score.fasta.part$P") {
					open (IN, "$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.score.fasta.part$P") or die;
					while (<IN>) {
						chomp;
						if ($_=~s/^>//) {
							my ($pos,$num_score) = (split "#",$_)[1,3];
							$pam = (split "#",$_)[2];
							my ($num,$score) = split "_",$num_score;
							print REI "Intergenic\t$pos\t$num\t$score\t";
						}else{
							print REI "$_"."$pam\n";
						}
					}
					close IN;
					unlink ("$outdir/$filename.User.sgRNA/$filename.$chr.intergenic.sgRNA.score.fasta.part$P") or die;
				}
			}
		}
		close REG;
		close REI;
	}

#############################################FASTQ FILE###############################################

	if ($input =~ /(fq.gz|fastq.gz)$/i) {
		print LOG  "# Your data $input is fastq.gz format.\n";
		print  "# Your data $input is fastq.gz format.\n";
		print LOG  "# Extracting possible sgRNA from fastq file.\n";
		print  "# Extracting possible sgRNA from fastq file.\n";
		my %hash;my $PAM;my $Reads;
		open(FASTQ,"gzip -dc $input|") or die;
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't open $filename.gene.sgRNA.fasta for writing!\n";
		while (<FASTQ>) {
			if(/^@/){
				chomp;
				my $position=tell(FASTQ);
				my $line_1=<FASTQ>;
				while ($line_1 =~ /(?=(\w{25}GG\w{3}))/g) {
					$Reads = $1;
					unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$hash{$Reads}++;
					}
				}
				$line_1 = reverse $line_1;
				$line_1 =~ tr/AGCT/TCGA/;
				while ($line_1 =~ /(?=(\w{25}GG\w{3}))/g) {
					$Reads = $1;
					unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$hash{$Reads}++;
					}
				}
				seek(FASTQ,$position, 0);
			}
		}
		close (FASTQ);
		foreach my $key (sort keys %hash) {
			$PAM = substr($key,24,3);
			print OUT ">$filename\_sgRNA#$PAM\#$hash{$key}\n$key\n";
		}
		close OUT;undef %hash;
		my $line = `wc -l $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta`;
		if ($line =~ /^(\d+)/){
			$line = $1;
			$line = $line / 2; 
		}
		my $file_line = (int($line/$process) + 1) * 2;
		system("split -d -l $file_line $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part");
		unlink("$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't delete $filename.sgRNA.fasta file";

		print LOG  "# Calculates the Rule set 2 score for the sgRNA.\n";
		print  "# Calculates the Rule set 2 score for the sgRNA.\n";
		for (my $P=0;$P<$process;$P++) {
			$pm->start and next;
			$P = sprintf("%02d",$P);
			my @fqgz_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part$P","--output","$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			system(@fqgz_rs2);
			if($? == -1) {
				die "system @fqgz_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		print LOG  "# Output the user's database.\n";
		print  "# Output the user's database.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.gene.sgRNA.db.fastq.txt");
		for (my $P=0;$P<$process;$P++) {
			$P = sprintf("%02d",$P);
			open (IN,"$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			my $pam;
			while (<IN>) {
				chomp;
				if ($_=~s/^>//) {
					my ($name,$num_score) = (split "#",$_)[0,2];
					$pam = (split "#",$_)[1];
					my ($num,$score) = split "_",$num_score;
					print OUT "$name\t$num\t$score\t";
				}else{
					print OUT "$_"."$pam\n";
				}
			}
			close IN;
			unlink ("$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P") or die;
		}
		close OUT;
	}elsif ($input =~ /(fq|fastq)$/i) {
		print LOG  "# Your data $input is fastq format.\n";
		print  "# Your data $input is fastq format.\n";
		print LOG  "# Extracting possible sgRNA from fastq file.\n";
		print  "# Extracting possible sgRNA from fastq file.\n";
		my %hash;my $PAM;my $Reads;
		open(FASTQ,"$input") or die;
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't open $filename.gene.sgRNA.fasta for writing!\n";
		while (<FASTQ>) {
			if(/^@/){
				chomp;
				my $position=tell(FASTQ);
				my $line_1=<FASTQ>;
				while ($line_1 =~ /(?=(\w{25}GG\w{3}))/g) {
					$Reads = $1;
					unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$hash{$Reads}++;
					}
				}
				$line_1 = reverse $line_1;
				$line_1 =~ tr/AGCT/TCGA/;
				while ($line_1 =~ /(?=(\w{25}GG\w{3}))/g) {
					$Reads = $1;
					unless ($Reads=~/N|T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$hash{$Reads}++;
					}
				}
				seek(FASTQ,$position, 0);
			}
		}
		close (FASTQ);
		foreach my $key (sort keys %hash) {
			$PAM = substr($key,24,3);
			print OUT ">$filename\_sgRNA#$PAM\#$hash{$key}\n$key\n";
		}
		close OUT;undef %hash;
		my $line = `wc -l $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta`;
		if ($line =~ /^(\d+)/){
			$line = $1;
			$line = $line / 2; 
		}
		my $file_line = (int($line/$process) + 1) * 2;
		system("split -d -l $file_line $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part");
		unlink("$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't delete $filename.sgRNA.fasta file";

		print LOG  "# Calculates the Rule set 2 score for the sgRNA.\n";
		print  "# Calculates the Rule set 2 score for the sgRNA.\n";
		for (my $P=0;$P<$process;$P++) {
			$pm->start and next;
			$P = sprintf("%02d",$P);
			my @fq_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part$P","--output","$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			system(@fq_rs2);
			if($? == -1) {
				die "system @fq_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		print LOG  "# Output the user's database.\n";
		print  "# Output the user's database.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.gene.sgRNA.db.fastq.txt");
		for (my $P=0;$P<$process;$P++) {
			$P = sprintf("%02d",$P);
			open (IN,"$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			my $pam;
			while (<IN>) {
				chomp;
				if ($_=~s/^>//) {
					my ($name,$num_score) = (split "#",$_)[0,2];
					$pam = (split "#",$_)[1];
					my ($num,$score) = split "_",$num_score;
					print OUT "$name\t$num\t$score\t";
				}else{
					print OUT "$_"."$pam\n";
				}
			}
			close IN;
			unlink ("$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P") or die;
		}
		close OUT;
	}	

#############################################FASTA FILE###############################################

	if ($input =~ /(fa.gz|fasta.gz)$/i) {
		print LOG  "# Your data $input is fasta.gz format.\n";
		print  "# Your data $input is fasta.gz format.\n";
		print LOG  "# Reading fasta file.\n";
		print  "# Reading fasta file.\n";
		my %contig;my($PAM,$i,$pos,$reads);
		open(FASTA,"gzip -dc $input|") or die;
		while (<FASTA>) {
			chomp;
			if ($_=~s/^>//) {
				$i=$_;
				$contig{$i}=undef;
			}else{
				$contig{$i}.=$_;
			}
		}
		close FASTA;
		print LOG  "# Extracting possible sgRNA from fasta file.\n";
		print  "# Extracting possible sgRNA from fasta file.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die;
		foreach my $contig_id (sort keys %contig) {
			my $length = length($contig{$contig_id});
			while ($contig{$contig_id}=~/(?=(\w{25}GG\w{3}))/g) {
				$pos = pos($contig{$contig_id}) + 4;
				$reads = $1;
				$PAM = substr($1,24,3);
				unless ($reads=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
					print OUT ">$contig_id\#+$pos\#$PAM\n$reads\n";
				}
			}
			$contig{$contig_id} = reverse $contig{$contig_id};
			$contig{$contig_id} =~ tr/ACGT/TGCA/;
			while ($contig{$contig_id}=~/(?=(\w{25}GG\w{3}))/g) {
				$pos = $length - pos($contig{$contig_id}) - 26;
				$reads = $1;
				$PAM = substr($1,24,3);
				unless ($reads=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
					print OUT ">$contig_id\#-$pos\#$PAM\n$reads\n";	
				}
			}
		}
		close OUT;
		my $line = `wc -l $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta`;
		if ($line =~ /^(\d+)/){
			$line = $1;
			$line = $line / 2; 
		}
		my $file_line = (int($line/$process) + 1) * 2;
		system("split -d -l $file_line $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part");
		unlink("$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't delete $filename.sgRNA.fasta file";

		print LOG  "# Calculates the Rule set 2 score for the sgRNA.\n";
		print  "# Calculates the Rule set 2 score for the sgRNA.\n";
		for (my $P=0;$P<$process;$P++) {
			$pm->start and next;
			$P = sprintf("%02d",$P);
			my @fagz_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part$P","--output","$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			system(@fagz_rs2);
			if($? == -1) {
				die "system @fagz_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		print LOG  "# Output the user's database.\n";
		print  "# Output the user's database.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.gene.sgRNA.db.fasta.txt");
		for (my $P=0;$P<$process;$P++) {
			$P = sprintf("%02d",$P);
			open (IN,"$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			my $pam;
			while (<IN>) {
				chomp;
				if ($_=~s/^>//) {
					my ($id,$pos,$PAM_score) = split "#",$_;
					$pam = (split "_",$PAM_score)[0];
					my ($score) = (split "_",$PAM_score)[1];
					print OUT "$id\t$pos\t$score\t";
				}else{
					print OUT "$_"."$pam\n";
				}
			}
			close IN;
			unlink ("$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P") or die;
		}
		close OUT;
	}elsif ($input =~ /(fa|fasta)$/i) {
		print LOG  "# Your data $input is fasta format.\n";
		print  "# Your data $input is fasta format.\n";
		print LOG  "# Reading fasta file.\n";
		print  "# Reading fasta file.\n";
		my %contig;my($PAM,$i,$pos,$reads);
		open(FASTA,"$input") or die;
		while (<FASTA>) {
			chomp;
			if ($_=~s/^>//) {
				$i=$_;
				$contig{$i}=undef;
			}else{
				$contig{$i}.=$_;
			}
		}
		close FASTA;
		print LOG  "# Extracting possible sgRNA from fasta file.\n";
		print  "# Extracting possible sgRNA from fasta file.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die;
		foreach my $contig_id (sort keys %contig) {
			my $length = length($contig{$contig_id});
			while ($contig{$contig_id}=~/(?=(\w{25}GG\w{3}))/g) {
				$pos = pos($contig{$contig_id}) + 4;
				$reads = $1;
				$PAM = substr($1,24,3);
				unless ($reads=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
					print OUT ">$contig_id\#+$pos\#$PAM\n$reads\n";
				}
			}
			$contig{$contig_id} = reverse $contig{$contig_id};
			$contig{$contig_id} =~ tr/ACGT/TGCA/;
			while ($contig{$contig_id}=~/(?=(\w{25}GG\w{3}))/g) {
				$pos = $length - pos($contig{$contig_id}) - 26;
				$reads = $1;
				$PAM = substr($1,24,3);
				unless ($reads=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
					print OUT ">$contig_id\#-$pos\#$PAM\n$reads\n";
				}
			}
		}
		close OUT;
		my $line = `wc -l $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta`;
		if ($line =~ /^(\d+)/){
			$line = $1;
			$line = $line / 2; 
		}
		my $file_line = (int($line/$process) + 1) * 2;
		system("split -d -l $file_line $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta $outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part");
		unlink("$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta") or die "Can't delete $filename.sgRNA.fasta file";

		print LOG  "# Calculates the Rule set 2 score for the sgRNA.\n";
		print  "# Calculates the Rule set 2 score for the sgRNA.\n";
		for (my $P=0;$P<$process;$P++) {
			$pm->start and next;
			$P = sprintf("%02d",$P);
			my @fa_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$outdir/$filename.User.sgRNA/$filename.sgRNA.fasta.part$P","--output","$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			system(@fa_rs2);
			if($? == -1) {
				die "system @fa_rs2 failed: $?";
			}
			$pm->finish;
		}
		$pm->wait_all_children;
		print LOG  "# Output the user's database.\n";
		print  "# Output the user's database.\n";
		open (OUT,">$outdir/$filename.User.sgRNA/$filename.gene.sgRNA.db.fasta.txt");
		for (my $P=0;$P<$process;$P++) {
			$P = sprintf("%02d",$P);
			open (IN,"$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P");
			my $pam;
			while (<IN>) {
				chomp;
				if ($_=~s/^>//) {
					my ($id,$pos,$PAM_score) = split "#",$_;
					$pam = (split "_",$PAM_score)[0];
					my ($score) = (split "_",$PAM_score)[1];
					print OUT "$id\t$pos\t$score\t";
				}else{
					print OUT "$_"."$pam\n";
				}
			}
			close IN;
			unlink ("$outdir/$filename.User.sgRNA/$filename.sgRNA.score.fasta.part$P") or die;
		}
		close OUT;
	}
print LOG "Your job is done.\n";
print LOG "################################# END ###########################################"."\n\n";
close LOG;
}
my $oo = time() - $time;
print "Total time consumption is $oo second.\n";

print "Your job is done.". "\n";
print  "################################# END ###########################################"."\n\n";
sub err {
	my ( $msg ) = @_;
	print "\n$msg\n\n";
	die "\n";
}
