#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
use Cwd;
use Parallel::ForkManager;


local $SIG{__WARN__} = sub {
	my $message = shift;
	die $message;
};
my $time = time();

my ($opt_mode,$Genome, $GFF3, $path, $label,$opt_help,$opt_up,$opt_down,$process,$spacer,$PAM_type);

GetOptions(	'm=s' => \$opt_mode,		#sgRNA designing mode(Cas9,Cpf1,custom)
			'i=s' => \$Genome,			#Input reference genome file (FASTA)
			'g=s' => \$GFF3,			#The reference annotation file (GFF3)
			'o=s' => \$path,			#Output path
			'l=s' => \$label,			#Name prefix for output file (default: CRISPR)
			'h!'  => \$opt_help,		#Help message
			'U:i' => \$opt_up,			#The number of bases to use for extending exon upstream sequence (default: 0)
			'D:i' => \$opt_down,		#The number of bases to use for extending exon downstream sequence (default: 0)
			'p:i' => \$process,			#The process to use
			't:s' => \$PAM_type,		#The type of PAM (default: 5'-TTTN)
			'x:i' => \$spacer,			#Length of protospacer (default: 24 nt)
          );

$opt_mode ||= "Cas9";
my $dir_default = getcwd;             #default output
$path ||= $dir_default;
$label ||= "CRISPR";
$process ||= "1";
my $dir =$path;
my $errmsg = "Use '$0 -h' for quick help; for more information, please see README.";

my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,

For Cas9 mode:

	perl $0 -m cas9 -i ZmB73_v4.fa -g ZmB73_v4.gff3 -o /your_dir/ -l ZmB73 -U 15 -D 3 -p 10

For Cpf1 mode:

	perl $0 -m cpf1 -i ZmB73_v4.fa -g ZmB73_v4.gff3 -o /your_dir/ -l ZmB73 -t TTTV -p 10 -x 23

For Custom mode:

	perl $0 -m custom -i ZmB73_v4.fa -g ZmB73_v4.gff3 -o /your_dir/ -l ZmB73 -t NRG -p 10 -x 20

--- Required ---
	
	-m <string>	:Given specific PAM for sgRNA design : Cas9, Cpf1 or Custom (default: Cas9)

	  Cas9		:On-target: 20 nt protospacer + NGG, off-target: 20 nt + NRG
	  Cpf1		:On-target: TTTN/TTN + 23/24/25 nt protospacer
	  Custom	:On-target: 15-25 nt protospacer + custom PAM sequence

	-i <string>	:reference genome sequence file
	-g <string>	:reference genome annotation file

--- Options ---

	-o <string>	:Output path (default: $dir_default)
	-l <string>	:Name prefix for output file (default: CRISPR)
	-p <int>	:The number of process to use (default:1)

	For Cas9:

	 -U <int>	:An integer specifying the number of base pairs expanding 5'-end of exon boundary (default: 0, >15 is not allowed)
	 -D <int>	:An integer specifying the number of base pairs expanding 3'-end of exon boundary (default: 0, >3 is not allowed)

	For Cpf1 or Custom:

	  -x <int>	:Cpf1: Length of spacer: between 23 to 25 (default: 24 nt);
			:Custom: Length of spacer: between 15 to 25 (default: 20 nt);
	  -t <string>	:Cpf1: Type of PAM sequence: TTX or TTTX (X: One of A C G T R Y M K S W H B V D N; default: TTTN)
			:Custom: Type of PAM sequence (default: NGG)

	-h: show this help
};

err( $helpmsg ) if $opt_help;

$opt_mode = lc $opt_mode;
err( "ERROR - Unknown mode: -m $opt_mode (supported mode: Cas9, Cpf1, Custom).\n$errmsg" ) unless grep /^$opt_mode$/, qw( cas9 cpf1 custom );

if ( not ($GFF3 && $Genome) ) {
	print
		"\n" .
		"Please input the option, and press <Enter>.\n" .
		"For quick help, input '-h' and press <Enter>.\n" .
		"\n";
	die "\n";
}

if ($opt_mode eq "cas9") {
	$opt_up ||= "0";
	$opt_down ||= "0";
	if ($opt_up) {
		err("ERROR - Wrong option of number of bases allowed: -U $opt_up (supported: 0-15).\n$errmsg") if ($opt_up < 0 || $opt_up >15);
	}
	if ($opt_down) {
		err ("ERROR - Wrong option of number of bases allowed: -D $opt_down (supported: 0-3).\n$errmsg") if ($opt_down < 0 || $opt_down >3);
	}
}elsif ($opt_mode eq "cpf1") {
	$PAM_type ||= "TTTN";
	$PAM_type = uc $PAM_type;
	$spacer ||= "24";
	err( "ERROR - Unknown PAM type: -t $PAM_type (supported PAM type: TTTX, TTX).\n$errmsg" ) 
		unless grep /^$PAM_type$/, qw( TTTA TTTC TTTG TTTT TTTR TTTY TTTM TTTK TTTS TTTW TTTH TTTB TTTV TTTD TTTN TTA TTC TTG TTT TTR TTY TTM TTK TTS TTW TTH TTB TTV TTD TTN );
	err( "ERROR - Wrong option of number of length of protospacer allowed: -x $spacer (supported: 23-25).\n$errmsg" ) if ($spacer < 23 || $spacer > 25);
}elsif ($opt_mode eq "custom") {
	$PAM_type ||= "NGG";
	$PAM_type = uc $PAM_type;
	$spacer ||= "20";
	err( "ERROR - Wrong option of number of length of protospacer allowed: -x $spacer (supported: 15-25).\n$errmsg" ) if ($spacer < 15 || $spacer > 25);
}

my $PAM_len;
if ($PAM_type) {
	err( "ERROR - Unknown bases in PAM sequence: -t $PAM_type (supported bases: A C G T R Y M K S W H B V D N).\n$errmsg") if ($PAM_type =~ /[^ACGTRYMKSWHBVDN]/);
	$PAM_len = length ($PAM_type);
	$PAM_type =~ s/R/[AG]/g;
	$PAM_type =~ s/Y/[CT]/g;
	$PAM_type =~ s/M/[AC]/g;
	$PAM_type =~ s/K/[GT]/g;
	$PAM_type =~ s/S/[GC]/g;
	$PAM_type =~ s/W/[AT]/g;
	$PAM_type =~ s/H/[ATC]/g;
	$PAM_type =~ s/B/[GTC]/g;
	$PAM_type =~ s/V/[GAC]/g;
	$PAM_type =~ s/D/[GAT]/g;
	$PAM_type =~ s/N/[ATCG]/g;
}

#In order for this step to work correctly, the chromosome names (sequence headers) in the FASTA reference genome file must be the same as the first column of GFF3 annotation file
#e.g., >1, >2, or >3 corresponding 1, 2, or 3.

system("mkdir -p -m 755 $dir/$label.Potential.off.target");
system("mkdir -p -m 755 $dir/$label.Possible.sgRNA");
system("mkdir -p -m 755 $dir/$label.Seqmap.result");

open  (LOG, ">$dir/$label.Log.txt") || die "Can't open $label.Log.txt for writing!" ."\n";

print  LOG "######################################### Log #########################################". "\n\n";
#print "Writing Log information.                                                                      " ."\n";
print  LOG "#                                     CRISPR-Local                                        " ."\n";
print  LOG "#  ---A local single-guide RNA (sgRNA) design tool for non-reference plant genomes." ."\n";          
print  LOG "#                                                                                         " ."\n";
print  LOG "#         contact:  Jiamin Sun, MaizeGo, Email: sunjm0824\@webmail.hzau.edu.cn                    \n\n";
######



print "\n  Welcome to CRISPR-Local\n";
print " ---A local single-guide RNA (sgRNA) design tool for non-reference plant genomes.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 1.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjm0824\@webmail.hzau.edu.cn"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today     : $local_time\n\n";
print "# Designing mode:$opt_mode\n\n";
print LOG "# Time, begin at $local_time."."\n";
print "# Program RD-build: CRISPR sgRNA design by using reference genome.\n\n";
print LOG "# Program RD-bulid: CRISPR sgRNA design by using reference genome.\n";
print LOG "# Usage: perl $0 -m $opt_mode -i $Genome -g $GFF3 -o $path -l $label -p $process\n\n";

################################### Reading Fasta ###################################
print  "# Start program RD-build........\n";
print  "# Step1: Reading fasta file.\n";
print LOG  "# Start program RD-build........\n";
print LOG  "# Step1: Reading fasta file.\n";

open (FASTA,$Genome) or die "Can't open $Genome for reading!\n";
my %chrom; my $i;
while (<FASTA>) {
	chomp;
	if ($_=~/^(>\S+)/) {
		$i = $1;
		$i =~ s/(>|chr)//ig;
		$chrom{$i}=undef;
	}else{
		my $line = uc $_;
		$chrom{$i}.=$line;
	}
}
close FASTA;
print  LOG "# Genome fasta parsed.\n";
print  "# Genome fasta parsed.\n";

if ($opt_mode eq "cas9") {
	################################### Reading GFF3 ###################################
	
	my (%gene_picked, %gene_picked_rev, %gene_START, %gene_END, %gene_CHR, %trans_gene);
	my (%exon_picked, %exon_picked_rev, %exon_CHR, %exon_START, %exon_END, %count_exon, %info_exon, %exon_length, %trans_exon);
	my %count_mRNA;
	my %TSS;my $gene_id;
	print  LOG "# Step2: Reading GFF3 file.\n";
	print  "# Step2: Reading GFF3 file.\n";
	open (GFF3,">$dir/$label.tmp.sorted.gff3") or die "Can't open $label.temp.sorted.gff3 for writing!\n";
	open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
	while (<GFF>) {
		next if ($_=~/^(\#|\n)/);
		my @line = (split /\t/,$_);
		if ($line[2] =~/gene$/) {
			my$ID="";
			if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
				$ID = $1;
				$ID =~ s/(NAME|id)\=//i;
				$ID =~ s/(gene|transcript|exon|CDS)://i;
			}else{
				next;
			}
			$line[8] = "$ID;\n";
			my $des = join("\t",@line[0..8]);
			print GFF3 "$des";
		}
	}
	close GFF;
	my @feature = qw/RNA transcript exon CDS/;
	foreach my $key (@feature) {
		open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
		while (<GFF>) {
			next if ($_=~/^(\#|\n)/);
			my @line = (split /\t/,$_);
			if ($line[2] =~/$key$/) {
				my$parent="";my$ID="";
				if ($line[8] =~ /(PARENT\=\S+?)(\n|;)/i) {
					$parent = $1;
					$parent =~ s/PARENT\=//i;
					$parent =~ s/(gene|transcript)://i;
				}else{
					next;
				}
				if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
					$ID = $1;
					$ID =~ s/(NAME|id)\=//i;
					$ID =~ s/(gene|transcript|exon|CDS)://i;
				}else{
					next;
				}
				$line[8] = "$ID;$parent\n";
				my $des = join("\t",@line[0..8]);
				print GFF3 "$des";
			}
		}
		close GFF;
	}
	close GFF3;
	open (GFF3,"$dir/$label.tmp.sorted.gff3") or die "Can't open $GFF3.sorted for writing!\n";
	open (INDEX,">$dir/$label.gene.index") or die;
	while (<GFF3>) {
		next if ($_=~/^\#/) ;
		chomp $_;
		my($chr,undef,$type,$start,$end,undef,$strand,undef,$attrs)= (split /\t/,$_);
		$chr =~s/chr//ig;
		my @attr = split( ";", $attrs );
		if (exists $chrom{$chr}) {
			if ($type =~/gene$/) {
				my $gene_name = $attr[0];
				my $gene_start = $start;
				my $gene_end = $end;
				print INDEX "$gene_name\t$chr:$strand"."$gene_start..$gene_end\n";
				my $gene_seq = substr($chrom{$chr}, $gene_start-1, $gene_end-$gene_start+1);
				$gene_START{$gene_name} = int($gene_start);
				$gene_END{$gene_name} = int($gene_end);
				$gene_CHR{$gene_name} = $chr;
				$gene_picked{$gene_name}=$gene_seq;
				$gene_seq = reverse $gene_seq;
				$gene_seq =~ tr/ACGT/TGCA/;
				$gene_picked_rev{$gene_name}=$gene_seq;
			}elsif ($type =~/exon/) {
				$gene_id = $trans_gene{$attr[1]};
				$info_exon{$gene_id}{$attr[0]}.= $attr[1].";";
				$count_exon{$gene_id}{$attr[0]}++;
				my $exon_name = $attr[0];
				my $exon_start = $start;
				my $exon_end = $end;
				$trans_exon{$attr[1]}{$exon_name} = "";
				my $exon_seq = substr($chrom{$chr}, $exon_start-$opt_up-1-4, int(($exon_end-$exon_start+$opt_up+$opt_down+1+4+3)));
				my $exon_seq_rev = substr($chrom{$chr}, $exon_start-$opt_down-1-3, int(($exon_end-$exon_start+$opt_up+$opt_down+1+4+3)));
				$exon_CHR{$exon_name} = $chr;
				$exon_START{$exon_name} = int($exon_start);
				$exon_END{$exon_name} = int($exon_end);
				$exon_length{$exon_name} = int(($exon_end-$exon_start+1));
				$exon_picked{$exon_name}=$exon_seq;
				$exon_seq_rev = reverse $exon_seq_rev;
				$exon_seq_rev =~ tr/ACGT/TGCA/;
				$exon_picked_rev{$exon_name}=$exon_seq_rev;
			}elsif($type =~ /RNA$|transcript/) {
				$count_mRNA{$attr[1]}++;
				$trans_gene{$attr[0]}=$attr[1];
			}elsif($type =~ /CDS$/) {
				if (not exists $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
				if ($start < $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
			}
		}
	}
	close GFF3;
	close INDEX;
	unlink("$dir/$label.tmp.sorted.gff3") or die "Can't delete tmp.sorted.gff3 file";
	print LOG  "# GFF3 file parsed.\n";
	print  "# GFF3 file parsed.\n";
	
	################################### Extract POT sites ###################################
	
	print LOG  "# Step3: Extracting potential off-target sites.\n";
	print  "# Step3: Extracting potential off-target sites.\n";
	
	my ($p, $OTseq, $pam, $pos);
	open (OUT,">$dir/$label.Potential.off.target/$label.POT.fasta") or die "Can't open $label.POT.fasta for writing!\n";
	foreach my $gene (sort keys %gene_picked) {
		while ($gene_picked{$gene}=~/(?=((\w{20}).(AG|GG)))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,0,20);
				$pam = substr($1,20,3);
				$pos = pos($gene_picked{$gene});
				$p = $gene_START{$gene} + $pos;
				print OUT ">$gene\#$gene_CHR{$gene}:+$p\#$pam\n$OTseq\n";
			}
		}
		$p=0;
		while ($gene_picked_rev{$gene}=~/(?=((\w{20}).(AG|GG)))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,0,20);
				$pam = substr($1,20,3);
				$pos = pos($gene_picked_rev{$gene});
				$p = $gene_END{$gene} - $pos -22;
				print OUT ">$gene\#$gene_CHR{$gene}:-$p\#$pam\n$OTseq\n";
			}
		}
		$p=0;
	}
	close OUT;
	print LOG  "# Potential off-target sites extracted.\n";
	print  "# Potential off-target sites extracted.\n";
	
	################################### Extract possible sgRNA ###################################
	
	print LOG  "# Step4: Extracting possible sgRNA.\n";
	print  "# Step4: Extracting possible sgRNA.\n";
	
	my $P=1;my $read_cnt=0;
	my ($relative_pos, $absolute_pos);
	my ($reads, $PAM, $string, $EXread, $score, $attr1s);
	open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die "Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
	foreach my $transcript (sort keys %trans_exon) {
		foreach $attr1s (sort keys %{$trans_exon{$transcript}}) {
			if (exists $exon_picked{$attr1s}) {
				while ($exon_picked{$attr1s}=~/(?=(\w{25}GG\w{3}))/g) {
					$reads = substr($1,4,20);
					$PAM = substr($1,24,3);
					$EXread = $1;
					$pos = pos($exon_picked{$attr1s}) + 4 - 4 - $opt_up;
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_START{$attr1s} + $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
					
						print OUT1 "$string\#$exon_CHR{$attr1s}:+"."$absolute_pos\#$PAM\n$EXread\n";
						if ($read_cnt>=50000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die"Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
						}
					}
				}
				while ($exon_picked_rev{$attr1s}=~/(?=(\w{25}GG\w{3}))/g) {
					$reads = substr($1,4,20);
					$PAM = substr($1,24,3);
					$EXread = $1;
					$pos = pos($exon_picked_rev{$attr1s}) + 25 - 3 - $opt_up;
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_END{$attr1s} - $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
						
						print OUT1 "$string\#$exon_CHR{$attr1s}:-"."$absolute_pos\#$PAM\n$EXread\n";
						if ($read_cnt>=50000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die"Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
						}
					}
				}
			}
			delete $exon_picked{$attr1s};
			delete $exon_picked_rev{$attr1s};
		}
	}
	close OUT1;
	my $pm = Parallel::ForkManager->new($process);
	my $file_num = $P;
	if($process == -1){
		$pm->set_max_procs($file_num);
	}elsif($process == 0){
		$pm->set_max_procs(0);
	}elsif($file_num > $process){
		$pm->set_max_procs($process);
	}else{
		$pm->set_max_procs($file_num);
	}
	print LOG  "# Possible sgRNA extracted.\n";
	print  "# Possible sgRNA extracted.\n";
	
	################################### INDEX file ###################################
	
	print LOG  "# Step5: Creating index file.\n";
	print  "# Step5: Creating index file.\n";
	
	open (OUT2,">$dir/$label.annotation.index") or die "Can't open $label.annotation.index for writing!\n";
	print OUT2 "Gene_id\tExon_id\tPosition\tNum_exon\tNum_transcript\tTranscript_id\n";
	foreach my $key1 (sort keys %count_mRNA) {
		foreach my $key2 (sort keys %{$info_exon{$key1}}) {
			print OUT2 "$key1\t$key2\t$gene_CHR{$key1}:$exon_START{$key2}..$exon_END{$key2}\t$count_exon{$key1}{$key2}\t$count_mRNA{$key1}\t$info_exon{$key1}{$key2}\n";
		}
	}
	close OUT2;
	print LOG  "# Index file created.\n";
	print  "# Index file created.\n";
	
	#=pod
	###################################### sgRNA on-target score ###################################
	print LOG  "# Step6: Calculates the Rule set 2 score for the given 30-mer sgRNA.\n";
	print  "# Step6: Calculates the Rule set 2 score for the given 30-mer sgRNA.\n";
	
	for (my $j=1; $j<=$file_num; $j++) {
		if (-e "$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta") {
			my $line1 = `wc -l $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.fasta`;
			my $line2 = `wc -l $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta`;
			if ($line1) {
				$line1 =$& if ($line1 =~ /^(\d+)/);
			}
			if ($line2) {
				$line2 =$& if ($line2 =~ /^(\d+)/);
			}
			if ($line1 == $line2) {
				system("rm $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.fasta");
				next;
			}
		}
		$pm->start and next;
		my @args_rs2 = ("python","Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py","--input","$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.fasta","--output","$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta");
		LABEL1:{
			system(@args_rs2);
			if($? == -1) {
				die "system @args_rs2 failed: $?";
				redo LABEL1;
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	###################################### seqmap ###################################
	print  "# Step7: Mapping Possible sgRNA to potential off-target sites.\n";
	print LOG  "# Step7: Mapping Possible sgRNA to potential off-target sites.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my $args_seqmap = ("./seqmap-1.0.12-linux-64 4 $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta $dir/$label.Potential.off.target/$label.POT.fasta $dir/$label.Seqmap.result/$label.seqmap_output.$j.txt /output_all_matches /skip_N /no_duplicate_probes /forward_strand /silent >> $dir/$label.Seqmap.result/$label.seqmap.LOG.txt");
		LABEL2:{
			system($args_seqmap);
			if($? == -1) {
				die "system $args_seqmap failed: $?";
			}
			if (-e "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt") {
				$pm->finish;
				next;
			}else{
				redo LABEL2;
			}
		}
	}
	$pm->wait_all_children;
	############################## format seqmap result ###################################
	print  "# Step8: Format seqmap result file.\n";
	print LOG  "# Step8: Format seqmap result file.\n";
	
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my @args_format = ("perl","sgRNA_reformatting.pl", "$opt_mode","$dir/$label.annotation.index", "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt", "$dir/$label.Seqmap.result/$label.result.$j.txt");
		LABEL3:{
			system(@args_format);
			if($? == -1) {
				die "system @args_format failed: $?";
				redo LABEL3;
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	############################## Calculate CFD score ###################################
	print  "# Step9: Calculates the Cutting Frequency Determination score.\n";
	print LOG  "# Step9: Calculates the Cutting Frequency Determination score.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my @args_CFD = ("python","cfd-score-calculator.py","--input","$dir/$label.Seqmap.result/$label.result.$j.txt","--output","$dir/$label.Seqmap.result/$label.CFD.$j.txt");
		system(@args_CFD);
		if($? == -1) {
			die "system @args_CFD failed: $?";
		}
		$pm->finish;
	}
	$pm->wait_all_children;
	############################## Merging and filtering ###################################
	print  "# Step10: Merging and filtering the CFD result.\n";
	print LOG  "# Step10: Merging and filtering the CFD result.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		system("cat $dir/$label.Seqmap.result/$label.CFD.$j.txt >> $dir/$label.Seqmap.result/$label.CFD.txt");
		unlink("$dir/$label.Seqmap.result/$label.CFD.$j.txt") or die "Can't delete $label.CFD.$j.txt file";
	}
	
	my %sgRNA;my %score;
	open (CFD,"$dir/$label.Seqmap.result/$label.CFD.txt") or die;
	while (<CFD>) {
		chomp;
		my($gene,$position,$rs2,$CFD) = (split /\t/,$_)[0,1,3,7];
		if (exists $sgRNA{"$gene:$position"}{$rs2}) {
			if ($CFD > $score{"$gene:$position"}{$rs2}) {
				$sgRNA{"$gene:$position"}{$rs2}="$_";
				$score{$position}{$rs2}=$CFD;
			}
		}else{
			$sgRNA{"$gene:$position"}{$rs2}="$_";
			$score{"$gene:$position"}{$rs2}=$CFD;
		}
	}
	close CFD;
	open (OUT3,">$dir/$label.Seqmap.result/$label.filt.result.txt") or die;
	
	foreach my $key1 (sort keys %sgRNA) {
		foreach my $key2 (sort keys %{$sgRNA{$key1}}) {
			print OUT3 "$sgRNA{$key1}{$key2}\n";
		}
	}
	close OUT3;
	unlink("$dir/$label.Seqmap.result/$label.CFD.txt") or die "Can't delete $label.CFD.txt file";
	
	system("mkdir -p -m 755 $dir/$label.sgRNA.database");
	system("sort -t\$'\t' -k1,1 -k4nr,4 $dir/$label.Seqmap.result/$label.filt.result.txt > $dir/$label.sgRNA.database/$label.reference.database.tmp");
	unlink("$dir/$label.Seqmap.result/$label.filt.result.txt") or die "Can't delete $label.filt.result.txt file";
	open (TMP, "$dir/$label.sgRNA.database/$label.reference.database.tmp") or die;
	open (RD, ">$dir/$label.sgRNA.database/$label.reference.database.txt") or die;
	while (<TMP>){
		chomp;
		my @TMP = (split /\t/,$_);
		my $tmp_guide = substr($TMP[2],0,20);
		my $tmp_pam = substr($TMP[2],20,3);
		my $tmp_ot_guide = substr($TMP[6],0,20);
		my $tmp_ot_pam = substr($TMP[6],20,3);
		$TMP[2] = $tmp_guide."+".$tmp_pam;
		$TMP[6] = $tmp_ot_guide."+".$tmp_ot_pam;
		my $TMP_line = join("\t",@TMP[0..10]);
		print RD "$TMP_line\n";
	}
	close TMP;
	close RD;
	unlink("$dir/$label.sgRNA.database/$label.reference.database.tmp") or die "$dir/$label.sgRNA.database/$label.reference.database.tmp";
}elsif ($opt_mode eq "cpf1") {
	################################### Reading GFF3 ###################################

	my (%gene_picked, %gene_picked_rev, %gene_START, %gene_END, %gene_CHR, %trans_gene);
	my (%exon_picked, %exon_picked_rev, %exon_CHR, %exon_START, %exon_END, %count_exon, %info_exon, %exon_length, %trans_exon);
	my %count_mRNA;
	my %TSS;my $gene_id;
	print  LOG "# Step2: Reading GFF3 file.\n";
	print  "# Step2: Reading GFF3 file.\n";
	open (GFF3,">$dir/$label.tmp.sorted.gff3") or die "Can't open $label.temp.sorted.gff3 for writing!\n";
	open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
	while (<GFF>) {
		next if ($_=~/^(\#|\n)/);
		my @line = (split /\t/,$_);
		if ($line[2] =~/gene$/) {
			my$ID="";
			if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
				$ID = $1;
				$ID =~ s/(NAME|id)\=//i;
				$ID =~ s/(gene|transcript|exon|CDS)://i;
			}else{
				next;
			}
			$line[8] = "$ID;\n";
			my $des = join("\t",@line[0..8]);
			print GFF3 "$des";
		}
	}
	close GFF;
	my @feature = qw/RNA transcript exon CDS/;
	foreach my $key (@feature) {
		open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
		while (<GFF>) {
			next if ($_=~/^(\#|\n)/);
			my @line = (split /\t/,$_);
			if ($line[2] =~/$key$/) {
				my$parent="";my$ID="";
				if ($line[8] =~ /(PARENT\=\S+?)(\n|;)/i) {
					$parent = $1;
					$parent =~ s/PARENT\=//i;
					$parent =~ s/(gene|transcript)://i;
				}else{
					next;
				}
				if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
					$ID = $1;
					$ID =~ s/(NAME|id)\=//i;
					$ID =~ s/(gene|transcript|exon|CDS)://i;
				}else{
					next;
				}
				$line[8] = "$ID;$parent\n";
				my $des = join("\t",@line[0..8]);
				print GFF3 "$des";
			}
		}
		close GFF;
	}
	close GFF3;
	open (GFF3,"$dir/$label.tmp.sorted.gff3") or die "Can't open $GFF3.sorted for writing!\n";
	open (INDEX,">$dir/$label.gene.index") or die;
	while (<GFF3>) {
		next if ($_=~/^\#/) ;
		chomp $_;
		my($chr,undef,$type,$start,$end,undef,$strand,undef,$attrs)= (split /\t/,$_);
		$chr =~s/chr//ig;
		my @attr = split( ";", $attrs );
		if (exists $chrom{$chr}) {
			if ($type =~/gene$/) {
				my $gene_name = $attr[0];
				my $gene_start = $start;
				my $gene_end = $end;
				print INDEX "$gene_name\t$chr:$strand"."$gene_start..$gene_end\n";
				my $gene_seq = substr($chrom{$chr}, $gene_start-1, $gene_end-$gene_start+1);
				$gene_START{$gene_name} = int($gene_start);
				$gene_END{$gene_name} = int($gene_end);
				$gene_CHR{$gene_name} = $chr;
				$gene_picked{$gene_name}=$gene_seq;
				$gene_seq = reverse $gene_seq;
				$gene_seq =~ tr/ACGT/TGCA/;
				$gene_picked_rev{$gene_name}=$gene_seq;
			}elsif ($type =~/exon/) {
				$gene_id = $trans_gene{$attr[1]};
				$info_exon{$gene_id}{$attr[0]}.= $attr[1].";";
				$count_exon{$gene_id}{$attr[0]}++;
				my $exon_name = $attr[0];
				my $exon_start = $start;
				my $exon_end = $end;
				$trans_exon{$attr[1]}{$exon_name} = "";
				my $exon_seq = substr($chrom{$chr}, $exon_start-1-4, int(($exon_end-$exon_start+1+4+4)));
				my $exon_seq_rev = $exon_seq;
				$exon_CHR{$exon_name} = $chr;
				$exon_START{$exon_name} = int($exon_start);
				$exon_END{$exon_name} = int($exon_end);
				$exon_length{$exon_name} = int(($exon_end-$exon_start+1));
				$exon_picked{$exon_name}=$exon_seq;
				$exon_seq_rev = reverse $exon_seq_rev;
				$exon_seq_rev =~ tr/ACGT/TGCA/;
				$exon_picked_rev{$exon_name}=$exon_seq_rev;
			}elsif($type =~ /RNA$|transcript/) {
				$count_mRNA{$attr[1]}++;
				$trans_gene{$attr[0]}=$attr[1];
			}elsif($type =~ /CDS$/) {
				if (not exists $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
				if ($start < $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
			}
		}
	}
	close GFF3;
	close INDEX;
	unlink("$dir/$label.tmp.sorted.gff3") or die "Can't delete tmp.sorted.gff3 file";
	print LOG  "# GFF3 file parsed.\n";
	print  "# GFF3 file parsed.\n";
	
	################################### Extract POT sites ###################################
	
	print LOG  "# Step3: Extracting potential off-target sites.\n";
	print  "# Step3: Extracting potential off-target sites.\n";
	
	my ($p, $OTseq, $pam, $pos);
	open (OUT,">$dir/$label.Potential.off.target/$label.POT.fasta") or die "Can't open $label.POT.fasta for writing!\n";
	foreach my $gene (sort keys %gene_picked) {
		while ($gene_picked{$gene}=~/(?=($PAM_type\w{$spacer}))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,$PAM_len,$spacer);
				$pam = substr($1,0,$PAM_len);
				$pos = pos($gene_picked{$gene});
				$p = $gene_START{$gene} + $pos;
				print OUT ">$gene\#$gene_CHR{$gene}:+$p\#$pam\n$OTseq\n";
			}
		}	
		$p=0;
		while ($gene_picked_rev{$gene}=~/(?=($PAM_type\w{$spacer}))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,$PAM_len,$spacer);
				$pam = substr($1,0,$PAM_len);
				$pos = pos($gene_picked_rev{$gene});
				$p = $gene_END{$gene} - $pos - $spacer - $PAM_len + 1;
				print OUT ">$gene\#$gene_CHR{$gene}:-$p\#$pam\n$OTseq\n";
			}
		}
		$p=0;
	}
	close OUT;
	print LOG  "# Potential off-target sites extracted.\n";
	print  "# Potential off-target sites extracted.\n";
	
	################################### Extract possible sgRNA ###################################
	
	print LOG  "# Step4: Extracting possible sgRNA.\n";
	print  "# Step4: Extracting possible sgRNA.\n";
	
	my $P=1;my $read_cnt=0;
	my ($relative_pos, $absolute_pos);
	my ($reads, $PAM, $string, $EXread, $score, $attr1s);
	open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.txt") or die "Can't open $label.Possible.sgRNA.$P.txt for writing!\n";
	print OUT1 "sgRNA_id\tsgRNA_seq\tchromatin_accessibility\n";
	foreach my $transcript (sort keys %trans_exon) {
		foreach $attr1s (sort keys %{$trans_exon{$transcript}}) {
			if (exists $exon_picked{$attr1s}) {
				while ($exon_picked{$attr1s}=~/(?=(\w{4}$PAM_type\w{$spacer}\w{4}))/g) {
					$reads = substr($1,$PAM_len+4,$spacer);
					$PAM = substr($1,4,$PAM_len);
					$EXread = $1;
					$pos = pos($exon_picked{$attr1s}) + 4 - 4;
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_START{$attr1s} + $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
					
						print OUT1 "$string\#$exon_CHR{$attr1s}:+$absolute_pos\#$PAM\t$EXread\t1\n";
						if ($read_cnt>=50000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.txt") or die"Can't open $label.Possible.sgRNA.$P.txt for writing!\n";
								print OUT1 "sgRNA_id\tsgRNA_seq\tchromatin_accessibility\n";
						}
					}
				}
				while ($exon_picked_rev{$attr1s}=~/(?=(\w{4}$PAM_type\w{$spacer}\w{4}))/g) {
					$reads = substr($1,$PAM_len+4,$spacer);
					$PAM = substr($1,4,$PAM_len);
					$EXread = $1;
					$pos = pos($exon_picked_rev{$attr1s}) + $spacer + $PAM_len - 1;
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_END{$attr1s} - $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
						
						print OUT1 "$string\#$exon_CHR{$attr1s}:-$absolute_pos\#$PAM\t$EXread\t1\n";
						if ($read_cnt>=50000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.txt") or die"Can't open $label.Possible.sgRNA.$P.txt for writing!\n";
								print OUT1 "sgRNA_id\tsgRNA_seq\tchromatin_accessibility\n";
						}
					}
				}
			}
			delete $exon_picked{$attr1s};
			delete $exon_picked_rev{$attr1s};
		}
	}
	close OUT1;
	my $pm = Parallel::ForkManager->new($process);
	my $file_num = $P;
	if($process == -1){
		$pm->set_max_procs($file_num);
	}elsif($process == 0){
		$pm->set_max_procs(0);
	}elsif($file_num > $process){
		$pm->set_max_procs($process);
	}else{
		$pm->set_max_procs($file_num);
	}
	print LOG  "# Possible sgRNA extracted.\n";
	print  "# Possible sgRNA extracted.\n";
	
	################################### INDEX file ###################################
	
	print LOG  "# Step5: Creating index file.\n";
	print  "# Step5: Creating index file.\n";
	
	open (OUT2,">$dir/$label.annotation.index") or die "Can't open $label.annotation.index for writing!\n";
	print OUT2 "Gene_id\tExon_id\tPosition\tNum_exon\tNum_transcript\tTranscript_id\n";
	foreach my $key1 (sort keys %count_mRNA) {
		foreach my $key2 (sort keys %{$info_exon{$key1}}) {
			print OUT2 "$key1\t$key2\t$gene_CHR{$key1}:$exon_START{$key2}..$exon_END{$key2}\t$count_exon{$key1}{$key2}\t$count_mRNA{$key1}\t$info_exon{$key1}{$key2}\n";
		}
	}
	close OUT2;
	print LOG  "# Index file created.\n";
	print  "# Index file created.\n";
	
	#=pod
	###################################### sgRNA on-target score ###################################
	print LOG  "# Step6: Calculates the Seq-deepCpf1 score for the given sgRNA.\n";
	print  "# Step6: Calculates the Seq-deepCpf1 score for the given sgRNA.\n";
	
	for (my $j=1; $j<=$file_num; $j++) {
		if (-e "$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.txt") {
			my $line1 = `wc -l $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.txt`;
			my $line2 = `wc -l $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.txt`;
			if ($line1) {
				$line1 =$& if ($line1 =~ /^(\d+)/);
			}
			if ($line2) {
				$line2 =$& if ($line2 =~ /^(\d+)/);
			}
			if ($line1 == $line2) {
				system("rm $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.txt");
				next;
			}
		}
		$pm->start and next;
		my @args_rs2 = ("python","DeepCpf1_Code/DeepCpf1.py","$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.txt","$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.txt");
		LABEL1:{
			system(@args_rs2);
			if($? == -1) {
				die "system @args_rs2 failed: $?";
				redo LABEL1;
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	###################################### seqmap ###################################
	print  "# Step7: Mapping Possible sgRNA to potential off-target sites.\n";
	print LOG  "# Step7: Mapping Possible sgRNA to potential off-target sites.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		open (TXT,"$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.txt") or die;
		open (FAS,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta") or die;
		readline TXT;
		while (<TXT>){
			chomp;
			my ($Deep_id,$Deep_seq,$Deep_score) = (split "\t",$_)[0,1,3];
			$Deep_seq = substr($Deep_seq,4+$PAM_len,$spacer);
			print FAS "$Deep_id#$Deep_score\n$Deep_seq\n";
		}
		close TXT;close FAS;
		my $args_seqmap = ("./seqmap-1.0.12-linux-64 4 $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.score.fasta $dir/$label.Potential.off.target/$label.POT.fasta $dir/$label.Seqmap.result/$label.seqmap_output.$j.txt /output_all_matches /skip_N /no_duplicate_probes /forward_strand /silent >> $dir/$label.Seqmap.result/$label.seqmap.LOG.txt");
		LABEL2:{
			system($args_seqmap);
			if($? == -1) {
				die "system $args_seqmap failed: $?";
			}
			if (-e "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt") {
				$pm->finish;
				next;
			}else{
				redo LABEL2;
			}
		}
	}
	$pm->wait_all_children;
	############################## format seqmap result ###################################
	print  "# Step8: Format seqmap result file.\n";
	print LOG  "# Step8: Format seqmap result file.\n";
	
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my @args_format = ("perl","sgRNA_reformatting.pl", "$opt_mode", "$dir/$label.annotation.index", "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt", "$dir/$label.Seqmap.result/$label.result.$j.txt");
		LABEL3:{
			system(@args_format);
			if($? == -1) {
				die "system @args_format failed: $?";
				redo LABEL3;
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	############################## Merging and filtering ###################################
	print  "# Step9: Merging and filtering the Seqmap result.\n";
	print LOG  "# Step9: Merging and filtering the Seqmap result.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		system("cat $dir/$label.Seqmap.result/$label.result.$j.txt >> $dir/$label.Seqmap.result/$label.result.txt");
		unlink("$dir/$label.Seqmap.result/$label.result.$j.txt") or die "Can't delete $label.result.$j.txt";
	}
	
	my %sgRNA;my %OT;my %info;my %OT_0;
	open (RES,"$dir/$label.Seqmap.result/$label.result.txt") or die;
	while (<RES>) {
		chomp;
		my($gene,$position,$sgRNA_seq,$Deep_score,$OT_gene,$OT_pos,$mismatch,$info1,$info2) = (split /\t/,$_)[0,1,2,3,4,5,7,8,9];
		my $sgRNA_seq_guide = substr($sgRNA_seq,$PAM_len,$spacer);
		my $sgRNA_seq_PAM = substr($sgRNA_seq,0,$PAM_len);
		$sgRNA_seq = $sgRNA_seq_PAM."+".$sgRNA_seq_guide;
		if (exists $sgRNA{"$gene:$position"}) {
			if ($position ne $OT_pos) {
				$OT{"$gene:$position"}{$mismatch}++;
				if ($mismatch == 0) {
					$OT_0{"$gene:$position"}.="$OT_gene($OT_pos);"
				}
			}
		}else{
			$OT{"$gene:$position"}{0} = 0;
			$OT{"$gene:$position"}{1} = 0;
			$OT{"$gene:$position"}{2} = 0;
			$OT{"$gene:$position"}{3} = 0;
			$OT{"$gene:$position"}{4} = 0;
			$sgRNA{"$gene:$position"}="$gene\t$position\t$sgRNA_seq\t$Deep_score";
			$info{"$gene:$position"}="$info1\t$info2";
		}
	}
	close RES;
	system("mkdir -p -m 755 $dir/$label.sgRNA.database");
	open (OUT3,">$dir/$label.sgRNA.database/$label.reference.database.txt") or die;
	
	foreach my $key (sort keys %sgRNA) {
		my $type;
		if ($OT{$key}{0} == 0 and $OT{$key}{1} == 0 and $OT{$key}{2} == 0 and $OT{$key}{3} == 0 and $OT{$key}{4} == 0) {
			$type = "NM";
		}elsif ($OT{$key}{0} == 1) {
			$type = "U0";
		}elsif ($OT{$key}{0} > 1) {
			$type = "R0";
		}elsif ($OT{$key}{1} == 1) {
			$type = "U1";
		}elsif ($OT{$key}{1} > 1) {
			$type = "R1";
		}elsif ($OT{$key}{2} == 1) {
			$type = "U2";
		}elsif ($OT{$key}{2} > 1) {
			$type = "R2";
		}elsif ($OT{$key}{3} == 1) {
			$type = "U3";
		}elsif ($OT{$key}{3} > 1) {
			$type = "R3";
		}elsif ($OT{$key}{4} == 1) {
			$type = "U4";
		}elsif ($OT{$key}{4} > 1) {
			$type = "R4";
		}
		my $OT_num = $OT{$key}{0}+$OT{$key}{1}+$OT{$key}{2}+$OT{$key}{3}+$OT{$key}{4};
		if (exists $OT_0{$key}) {
			print OUT3 "$sgRNA{$key}\t$OT_num\t$type\t$OT{$key}{0},$OT{$key}{1},$OT{$key}{2},$OT{$key}{3},$OT{$key}{4}\t$OT_0{$key}\t$info{$key}\tNA\n";
		}else{
			print OUT3 "$sgRNA{$key}\t$OT_num\t$type\t$OT{$key}{0},$OT{$key}{1},$OT{$key}{2},$OT{$key}{3},$OT{$key}{4}\tNA\t$info{$key}\tNA\n";
		}
	}
	close OUT3;
}elsif ($opt_mode eq "custom") {
	################################### Reading GFF3 ###################################

	my (%gene_picked, %gene_picked_rev, %gene_START, %gene_END, %gene_CHR, %trans_gene);
	my (%exon_picked, %exon_picked_rev, %exon_CHR, %exon_START, %exon_END, %count_exon, %info_exon, %exon_length, %trans_exon);
	my %count_mRNA;
	my %TSS;my $gene_id;
	print  LOG "# Step2: Reading GFF3 file.\n";
	print  "# Step2: Reading GFF3 file.\n";
	open (GFF3,">$dir/$label.tmp.sorted.gff3") or die "Can't open $label.temp.sorted.gff3 for writing!\n";
	open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
	while (<GFF>) {
		next if ($_=~/^(\#|\n)/);
		my @line = (split /\t/,$_);
		if ($line[2] =~/gene$/) {
			my$ID="";
			if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
				$ID = $1;
				$ID =~ s/(NAME|id)\=//i;
				$ID =~ s/(gene|transcript|exon|CDS)://i;
			}else{
				next;
			}
			$line[8] = "$ID;\n";
			my $des = join("\t",@line[0..8]);
			print GFF3 "$des";
		}
	}
	close GFF;
	my @feature = qw/RNA transcript exon CDS/;
	foreach my $key (@feature) {
		open (GFF,$GFF3) or die "Can't open $GFF3 for reading!\n";
		while (<GFF>) {
			next if ($_=~/^(\#|\n)/);
			my @line = (split /\t/,$_);
			if ($line[2] =~/$key$/) {
				my$parent="";my$ID="";
				if ($line[8] =~ /(PARENT\=\S+?)(\n|;)/i) {
					$parent = $1;
					$parent =~ s/PARENT\=//i;
					$parent =~ s/(gene|transcript)://i;
				}else{
					next;
				}
				if ($line[8] =~ /((NAME|ID)\=\S+?)(\n|;)/i) {
					$ID = $1;
					$ID =~ s/(NAME|id)\=//i;
					$ID =~ s/(gene|transcript|exon|CDS)://i;
				}else{
					next;
				}
				$line[8] = "$ID;$parent\n";
				my $des = join("\t",@line[0..8]);
				print GFF3 "$des";
			}
		}
		close GFF;
	}
	close GFF3;
	open (GFF3,"$dir/$label.tmp.sorted.gff3") or die "Can't open $GFF3.sorted for writing!\n";
	open (INDEX,">$dir/$label.gene.index") or die;
	while (<GFF3>) {
		next if ($_=~/^\#/) ;
		chomp $_;
		my($chr,undef,$type,$start,$end,undef,$strand,undef,$attrs)= (split /\t/,$_);
		$chr =~s/chr//ig;
		my @attr = split( ";", $attrs );
		if (exists $chrom{$chr}) {
			if ($type =~/gene$/) {
				my $gene_name = $attr[0];
				my $gene_start = $start;
				my $gene_end = $end;
				print INDEX "$gene_name\t$chr:$strand"."$gene_start..$gene_end\n";
				my $gene_seq = substr($chrom{$chr}, $gene_start-1, $gene_end-$gene_start+1);
				$gene_START{$gene_name} = int($gene_start);
				$gene_END{$gene_name} = int($gene_end);
				$gene_CHR{$gene_name} = $chr;
				$gene_picked{$gene_name}=$gene_seq;
				$gene_seq = reverse $gene_seq;
				$gene_seq =~ tr/ACGT/TGCA/;
				$gene_picked_rev{$gene_name}=$gene_seq;
			}elsif ($type =~/exon/) {
				$gene_id = $trans_gene{$attr[1]};
				$info_exon{$gene_id}{$attr[0]}.= $attr[1].";";
				$count_exon{$gene_id}{$attr[0]}++;
				my $exon_name = $attr[0];
				my $exon_start = $start;
				my $exon_end = $end;
				$trans_exon{$attr[1]}{$exon_name} = "";
				my $exon_seq = substr($chrom{$chr}, $exon_start-1, int(($exon_end-$exon_start+1)));
				my $exon_seq_rev = $exon_seq;
				$exon_CHR{$exon_name} = $chr;
				$exon_START{$exon_name} = int($exon_start);
				$exon_END{$exon_name} = int($exon_end);
				$exon_length{$exon_name} = int(($exon_end-$exon_start+1));
				$exon_picked{$exon_name}=$exon_seq;
				$exon_seq_rev = reverse $exon_seq_rev;
				$exon_seq_rev =~ tr/ACGT/TGCA/;
				$exon_picked_rev{$exon_name}=$exon_seq_rev;
			}elsif($type =~ /RNA$|transcript/) {
				$count_mRNA{$attr[1]}++;
				$trans_gene{$attr[0]}=$attr[1];
			}elsif($type =~ /CDS$/) {
				if (not exists $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
				if ($start < $TSS{$attr[1]}) {
					$TSS{$attr[1]}=$start;
				}
			}
		}
	}
	close GFF3;
	close INDEX;
	unlink("$dir/$label.tmp.sorted.gff3") or die "Can't delete tmp.sorted.gff3 file";
	print LOG  "# GFF3 file parsed.\n";
	print  "# GFF3 file parsed.\n";
	
	################################### Extract POT sites ###################################
	
	print LOG  "# Step3: Extracting potential off-target sites.\n";
	print  "# Step3: Extracting potential off-target sites.\n";
	
	my ($p, $OTseq, $pam, $pos);
	open (OUT,">$dir/$label.Potential.off.target/$label.POT.fasta") or die "Can't open $label.POT.fasta for writing!\n";
	foreach my $gene (sort keys %gene_picked) {
		while ($gene_picked{$gene}=~/(?=(\w{$spacer}$PAM_type))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,0,$spacer);
				$pam = substr($1,$spacer,$PAM_len);
				$pos = pos($gene_picked{$gene});
				$p = $gene_START{$gene} + $pos;
				print OUT ">$gene\#$gene_CHR{$gene}:+$p\#$pam\n$OTseq\n";
			}
		}	
		$p=0;
		while ($gene_picked_rev{$gene}=~/(?=(\w{$spacer}$PAM_type))/g) {
			unless ($1=~/[^AGCT]/) {
				$OTseq = substr($1,0,$spacer);
				$pam = substr($1,$spacer,$PAM_len);
				$pos = pos($gene_picked_rev{$gene});
				$p = $gene_END{$gene} - $pos - $spacer - $PAM_len + 1;
				print OUT ">$gene\#$gene_CHR{$gene}:-$p\#$pam\n$OTseq\n";
			}
		}
		$p=0;
	}
	close OUT;
	print LOG  "# Potential off-target sites extracted.\n";
	print  "# Potential off-target sites extracted.\n";
	
	################################### Extract possible sgRNA ###################################
	
	print LOG  "# Step4: Extracting possible sgRNA.\n";
	print  "# Step4: Extracting possible sgRNA.\n";
	
	my $P=1;my $read_cnt=0;
	my ($relative_pos, $absolute_pos);
	my ($reads, $PAM, $string, $EXread, $score, $attr1s);
	open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die "Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
	foreach my $transcript (sort keys %trans_exon) {
		foreach $attr1s (sort keys %{$trans_exon{$transcript}}) {
			if (exists $exon_picked{$attr1s}) {
				while ($exon_picked{$attr1s}=~/(?=(\w{$spacer}$PAM_type))/g) {
					$reads = substr($1,0,$spacer);
					$PAM = substr($1,$spacer,$PAM_len);
					$EXread = $1;
					$pos = pos($exon_picked{$attr1s});
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_START{$attr1s} + $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
					
						print OUT1 "$string\#$exon_CHR{$attr1s}:+$absolute_pos\#$PAM\#NA\n$reads\n";
						if ($read_cnt>=100000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die"Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
						}
					}
				}
				while ($exon_picked_rev{$attr1s}=~/(?=(\w{$spacer}$PAM_type))/g) {
					$reads = substr($1,0,$spacer);
					$PAM = substr($1,$spacer,$PAM_len);
					$EXread = $1;
					$pos = pos($exon_picked_rev{$attr1s}) + $spacer + $PAM_len - 1;
					unless ($EXread=~/[^AGCT]/ or $reads=~/T{4,20}|A{5,20}|C{5,20}|G{5,20}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|(GC){6,10}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/) {
						$read_cnt++;
						$string = ">$attr1s";
						$absolute_pos = int($exon_END{$attr1s} - $pos);
						if (exists $TSS{$transcript}) {
							$relative_pos = int($pos + $exon_START{$attr1s} - $TSS{$transcript});
							$string = "$string\#[$TSS{$transcript}:$exon_START{$attr1s}:$exon_length{$attr1s}:$pos:$relative_pos]";
						}else{
							$string = "$string\#[$exon_START{$attr1s}:$exon_length{$attr1s}:$pos]";
						}
						
						print OUT1 "$string\#$exon_CHR{$attr1s}:-$absolute_pos\#$PAM\#NA\n$reads\n";
						if ($read_cnt>=50000) {
							$read_cnt=0;
							close OUT1;
							$P++;
							open (OUT1,">$dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$P.fasta") or die"Can't open $label.Possible.sgRNA.$P.fasta for writing!\n";
						}
					}
				}
			}
			delete $exon_picked{$attr1s};
			delete $exon_picked_rev{$attr1s};
		}
	}
	close OUT1;
	my $pm = Parallel::ForkManager->new($process);
	my $file_num = $P;
	if($process == -1){
		$pm->set_max_procs($file_num);
	}elsif($process == 0){
		$pm->set_max_procs(0);
	}elsif($file_num > $process){
		$pm->set_max_procs($process);
	}else{
		$pm->set_max_procs($file_num);
	}
	print LOG  "# Possible sgRNA extracted.\n";
	print  "# Possible sgRNA extracted.\n";
	
	################################### INDEX file ###################################
	
	print LOG  "# Step5: Creating index file.\n";
	print  "# Step5: Creating index file.\n";
	
	open (OUT2,">$dir/$label.annotation.index") or die "Can't open $label.annotation.index for writing!\n";
	print OUT2 "Gene_id\tExon_id\tPosition\tNum_exon\tNum_transcript\tTranscript_id\n";
	foreach my $key1 (sort keys %count_mRNA) {
		foreach my $key2 (sort keys %{$info_exon{$key1}}) {
			print OUT2 "$key1\t$key2\t$gene_CHR{$key1}:$exon_START{$key2}..$exon_END{$key2}\t$count_exon{$key1}{$key2}\t$count_mRNA{$key1}\t$info_exon{$key1}{$key2}\n";
		}
	}
	close OUT2;
	print LOG  "# Index file created.\n";
	print  "# Index file created.\n";
	
	#=pod
	
	###################################### seqmap ###################################
	print  "# Step6: Mapping Possible sgRNA to potential off-target sites.\n";
	print LOG  "# Step6: Mapping Possible sgRNA to potential off-target sites.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my $args_seqmap = ("./seqmap-1.0.12-linux-64 4 $dir/$label.Possible.sgRNA/$label.Possible.sgRNA.$j.fasta $dir/$label.Potential.off.target/$label.POT.fasta $dir/$label.Seqmap.result/$label.seqmap_output.$j.txt /output_all_matches /skip_N /no_duplicate_probes /forward_strand /silent >> $dir/$label.Seqmap.result/$label.seqmap.LOG.txt");
		LABEL2:{
			system($args_seqmap);
			if($? == -1) {
				die "system $args_seqmap failed: $?";
			}
			if (-e "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt") {
				$pm->finish;
				next;
			}else{
				redo LABEL2;
			}
		}
	}
	$pm->wait_all_children;
	############################## format seqmap result ###################################
	print  "# Step7: Format seqmap result file.\n";
	print LOG  "# Step7: Format seqmap result file.\n";
	
	for (my $j=1; $j<=$file_num; $j++) {
		$pm->start and next;
		my @args_format = ("perl","sgRNA_reformatting.pl", "$opt_mode", "$dir/$label.annotation.index", "$dir/$label.Seqmap.result/$label.seqmap_output.$j.txt", "$dir/$label.Seqmap.result/$label.result.$j.txt");
		LABEL3:{
			system(@args_format);
			if($? == -1) {
				die "system @args_format failed: $?";
				redo LABEL3;
			}
			$pm->finish;
		}
	}
	$pm->wait_all_children;
	
	############################## Merging and filtering ###################################
	print  "# Step8: Merging and filtering the Seqmap result.\n";
	print LOG  "# Step8: Merging and filtering the Seqmap result.\n";
	for (my $j=1; $j<=$file_num; $j++) {
		system("cat $dir/$label.Seqmap.result/$label.result.$j.txt >> $dir/$label.Seqmap.result/$label.result.txt");
		unlink("$dir/$label.Seqmap.result/$label.result.$j.txt") or die "Can't delete $label.result.$j.txt";
	}
	
	my %sgRNA;my %OT;my %info;my %OT_0;
	open (RES,"$dir/$label.Seqmap.result/$label.result.txt") or die;
	while (<RES>) {
		chomp;
		my($gene,$position,$sgRNA_seq,$OT_gene,$OT_pos,$mismatch,$info1,$info2) = (split /\t/,$_)[0,1,2,4,5,7,8,9];
		my $sgRNA_seq_guide = substr($sgRNA_seq,0,$spacer);
		my $sgRNA_seq_PAM = substr($sgRNA_seq,$spacer,$PAM_len);
		$sgRNA_seq = $sgRNA_seq_guide."+".$sgRNA_seq_PAM;

		if (exists $sgRNA{"$gene:$position"}) {
			if ($position ne $OT_pos) {
				$OT{"$gene:$position"}{$mismatch}++;
				if ($mismatch == 0) {
					$OT_0{"$gene:$position"}.="$OT_gene($OT_pos);"
				}
			}
		}else{
			$OT{"$gene:$position"}{0} = 0;
			$OT{"$gene:$position"}{1} = 0;
			$OT{"$gene:$position"}{2} = 0;
			$OT{"$gene:$position"}{3} = 0;
			$OT{"$gene:$position"}{4} = 0;
			$sgRNA{"$gene:$position"}="$gene\t$position\t$sgRNA_seq\tNA";
			$info{"$gene:$position"}="$info1\t$info2";
		}
	}
	close RES;
	system("mkdir -p -m 755 $dir/$label.sgRNA.database");
	open (OUT3,">$dir/$label.sgRNA.database/$label.reference.database.txt") or die;
	
	foreach my $key (sort keys %sgRNA) {
		my $type;
		if ($OT{$key}{0} == 0 and $OT{$key}{1} == 0 and $OT{$key}{2} == 0 and $OT{$key}{3} == 0 and $OT{$key}{4} == 0) {
			$type = "NM";
		}elsif ($OT{$key}{0} == 1) {
			$type = "U0";
		}elsif ($OT{$key}{0} > 1) {
			$type = "R0";
		}elsif ($OT{$key}{1} == 1) {
			$type = "U1";
		}elsif ($OT{$key}{1} > 1) {
			$type = "R1";
		}elsif ($OT{$key}{2} == 1) {
			$type = "U2";
		}elsif ($OT{$key}{2} > 1) {
			$type = "R2";
		}elsif ($OT{$key}{3} == 1) {
			$type = "U3";
		}elsif ($OT{$key}{3} > 1) {
			$type = "R3";
		}elsif ($OT{$key}{4} == 1) {
			$type = "U4";
		}elsif ($OT{$key}{4} > 1) {
			$type = "R4";
		}
		my $OT_num = $OT{$key}{0}+$OT{$key}{1}+$OT{$key}{2}+$OT{$key}{3}+$OT{$key}{4};
		if (exists $OT_0{$key}) {
			print OUT3 "$sgRNA{$key}\t$OT_num\t$type\t$OT{$key}{0},$OT{$key}{1},$OT{$key}{2},$OT{$key}{3},$OT{$key}{4}\t$OT_0{$key}\t$info{$key}\tNA\n";
		}else{
			print OUT3 "$sgRNA{$key}\t$OT_num\t$type\t$OT{$key}{0},$OT{$key}{1},$OT{$key}{2},$OT{$key}{3},$OT{$key}{4}\tNA\t$info{$key}\tNA\n";
		}
	}
	close OUT3;
}
system("rm -rf $dir/$label.Seqmap.result");
system("rm -rf $dir/$label.Potential.off.target");
system("rm -rf $dir/$label.Possible.sgRNA");
#unlink("$dir/$label.annotation.index");

my $oo = time() - $time;
print "# Total time consumption is $oo second.\nDone!\n";

print "# Your job is done, open $label.reference.database.txt to check result.". "\n";
print  LOG "# Your job is done, open $label.reference.database.txt to check result.". "\n";

print  LOG "################################# END ###########################################"."\n\n";
close LOG;
#=cut

sub err {
	my ( $msg ) = @_;
	print "\n$msg\n\n";
	die "\n";
}
