#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
use Cwd;

local $SIG{__WARN__} = sub {
	my $message = shift;
	die $message;
};

my ($List, $RefDB, $UserDB, $Path, $opt_help,$label);
my $time = time();

GetOptions(	
			'g=s' => \$List,
			'i=s' => \$RefDB,
			'u:s' => \$UserDB,
			'o=s' => \$Path,
			'h!'  => \$opt_help,		#Help message
			'l:s' => \$label,
          );

my $dir_default = getcwd;             #default output
$Path ||= $dir_default;
my $dir =$Path;
my $errmsg = "Use '$0 -h' for quick help; for more information, please see README.";
$label ||= "DB_search";
my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,
	perl $0 -g ZmB73_query_gene.list -i ZmB73.reference.database.txt -o /your_dir/
	
	or

	perl $0 -g ZmB73_query_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/

--- Options ---

	-g <string>	:Query gene list file
	-i <string>	:Reference sgRNA database (RD)
	-u <string>	:User's sgRNA database (UD)
	-o <string>	:Output path (defualt: $dir_default)
	-l <label>	:Name prefix for output file (default:DB_search)

	-h		:show this help
};
err( $helpmsg ) if $opt_help;

err( "ERROR - No query gene list file provided.\n$errmsg" ) if !$List;
err( "ERROR - No reference database file provided.\n$errmsg" ) if !$RefDB;

print "\n  Welcome to CRISPR-Local\n";
print " Program DB_search: a local tool for search best target of query gene.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 3.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjm0824\@webmail.hzau.edu.cn"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today : $local_time\n\n";
print "# Program DB-search: Search the reference and user's database to select the sgRNA with high on-target score, low off-target score.\n";

print "# Reading gene list file.\n";
open (LIST,$List) or die "Can't open $List for reading!\n";
my %gene_list;
while (<LIST>) {
	chomp;
	$gene_list{$_}=undef;
}
my @list = keys %gene_list;
my $list = @list;
print "# There are $list gene in your list.\n";

my(%RD_sgRNA,%UD_sgRNA);
print "# Reading reference database file.\n";
open (SGRNA,$RefDB) or die $!;
while (<SGRNA>) {
	chomp;
	my($gene,$pos,$seq,$score,$OT_1,$OT_2,$OT_3,$OT_4,$OT_score)=(split /\t/,$_)[0,1,2,3,4,5,6,7,-1];
	if (exists $gene_list{$gene}) {
		$RD_sgRNA{$gene}{$seq}.="$gene\t$pos\t$seq\t$score\t$OT_1\t$OT_2\t$OT_3\t$OT_4\t$OT_score\n";
	}
}
close SGRNA;
my @RD_valid = sort keys %RD_sgRNA;
my $RD_num = @RD_valid;
my %gene_left = %gene_list;
print "# There are $RD_num available gene in the reference database.\n";
if ($list > $RD_num) {
	foreach my $key (@RD_valid) {
		delete $gene_left{$key};
	}
	open (LEFT,">$dir/$label\_Invalid_gene_RD.list") or die "Can't open $label\_Invalid_gene_RD.list for writing!\n";
	foreach my $key (sort keys %gene_left) {
		print LEFT "$key\n";
	}
	print "# The invalid gene is output to $label\_Invalid_gene_RD.list, please add annotation information of these invalid gene.\n";
	close LEFT;
}

if (!$UserDB) {
	print "# Output reference database search result.\n";
	open (OUT,">$dir/$label\_result_RD.txt") or die "Can't open $label\_result_RD.txt for writing!\n";
	foreach my $RD_valid (@RD_valid) {
		my @rs2 = sort keys %{$RD_sgRNA{$RD_valid}};
		for (my $i=0;$i<=$#rs2;$i++) {
			print OUT "$RD_sgRNA{$RD_valid}{$rs2[$i]}";
		}
	}
	close OUT;
	print "# Your job is done, check result.\n# END\n\n";
	my $oo = time() - $time;
	print "Total time consumption is $oo second.\nDone!\n";
	exit(0);
}

if ($UserDB =~ /gene.sgRNA.db.alignment.txt$/) {
	print "# Reading User's alignment database file.\n";
	open (UDB,$UserDB) or die;
	while (<UDB>) {
		chomp;
		my($gene,$pos,$num,$score,$seq) = split /\t/,$_;
		if (exists $gene_list{$gene}) {
			$UD_sgRNA{$gene}{$seq}.="$_\n";
		}
	}
	close UDB;
	my @UD_valid = sort keys %UD_sgRNA;
	my $UD_num = @UD_valid;
	%gene_left = %gene_list;
	print "# There are $UD_num available gene in the user's alignment database.\n";
	if ($list > $UD_num) {
		foreach my $key (@UD_valid) {
			delete $gene_left{$key};
		}
		open (LEFT,">$dir/$label\_Invalid_gene_UD.list") or die "Can't open $label\_Invalid_gene_UD.list for writing!\n";
		foreach my $key (sort keys %gene_left) {
			print LEFT "$key\n";
		}
		close LEFT;
		print "# The invalid gene is output to $label\_Invalid_gene_UD.list.\n";
	}
	print "# Find out the 'RO', 'UO' and 'BO' sgRNA.\n";
	my @inter = grep {$RD_sgRNA{$_}} @UD_valid;
	my %merge = map {$_ => 1} @RD_valid,@UD_valid;
	my @merge = sort keys (%merge);
	my @RDC_gene = grep {!$RD_sgRNA{$_}} @merge;
	my @UDC_gene = grep {!$UD_sgRNA{$_}} @merge;
	open (REF,">$dir/$label\_result_RO.txt");
	open (USER,">$dir/$label\_result_UO.txt");
	open (CORNA,">$dir/$label\_result_BO.txt");
	print "# Output the 'RO', 'UO' and 'BO' sgRNA result.\n";
	foreach my $RDC_gene (@RDC_gene) {
		my @rs2 = sort keys %{$UD_sgRNA{$RDC_gene}};
		for (my $i=0;$i<=$#rs2;$i++) {
			print USER "$UD_sgRNA{$RDC_gene}{$rs2[$i]}";
		}
	}
	foreach my $UDC_gene (@UDC_gene) {
		my @rs2 = sort keys %{$RD_sgRNA{$UDC_gene}};
		for (my $i=0;$i<=$#rs2;$i++) {
			print REF "$RD_sgRNA{$UDC_gene}{$rs2[$i]}";
		}
	}
	foreach my $inter (sort @inter) {
		my @RD_seq = sort keys %{$RD_sgRNA{$inter}};
		my @UD_seq = sort keys %{$UD_sgRNA{$inter}};
		my @inter_seq = grep {$RD_sgRNA{$inter}{$_}} @UD_seq;
		my %merge_seq = map {$_ => 1} @RD_seq,@UD_seq;
		my @merge_seq = sort keys (%merge_seq);
		my @RDC_seq = grep {!$RD_sgRNA{$inter}{$_}} @merge_seq;
		my @UDC_seq = grep {!$UD_sgRNA{$inter}{$_}} @merge_seq;

		for (my $i=0;$i<=$#RDC_seq;$i++) {
			print USER "$UD_sgRNA{$inter}{$RDC_seq[$i]}";
		}
		for (my $i=0;$i<=$#UDC_seq;$i++) {
			print REF "$RD_sgRNA{$inter}{$UDC_seq[$i]}";
		}
		for (my $i=0;$i<=$#inter_seq;$i++) {
			print CORNA "$RD_sgRNA{$inter}{$inter_seq[$i]}";
		}
	}
	close REF;close USER;close CORNA;
}elsif ($UserDB =~ /intergenic.sgRNA.db.alignment.txt$/) {
	print "# Reading User's intergenic alignment database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/$label\_intergenic_result_BO.txt");
	open (REF,">$dir/$label\_intergenic_result_RO.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_sgRNA{$RD_valid}}){
			$seq_num{$seq}=0;
		}
	}
	while (<UDB>) {
		chomp;
		my ($num,$inter_seq) = (split /\t/,$_)[2,-1];
		if (exists $seq_num{$inter_seq}) {
			$seq_num{$inter_seq}+=$num;
		}
	}
	close UDB;
	print "# Find out the 'RO' and 'BO' sgRNA.\n";
	foreach my $RD_valid (@RD_valid) {
		foreach my $key (keys %{$RD_sgRNA{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				print INTER "$RD_sgRNA{$RD_valid}{$key}";
			}else{
				print REF "$RD_sgRNA{$RD_valid}{$key}";
			}
		}
	}
	close INTER;close REF;
}elsif ($UserDB =~ /gene.sgRNA.db.fastq.txt$/) {
	print "# Reading User's fastq reads database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/$label\_result_BO.txt");
	open (REF,">$dir/$label\_result_RO.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_sgRNA{$RD_valid}}){
			$seq_num{$seq}=0;
		}
	}
	while (<UDB>) {
		chomp;
		my ($num,$inter_seq) = (split /\t/,$_)[1,-1];
		if (exists $seq_num{$inter_seq}) {
			$seq_num{$inter_seq}+=$num;
		}
	}
	close UDB;
	print "# Find out the 'RO' and 'BO' sgRNA.\n";
	foreach my $RD_valid (@RD_valid) {
		foreach my $key (keys %{$RD_sgRNA{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				print INTER "$RD_sgRNA{$RD_valid}{$key}";
			}else{
				print REF "$RD_sgRNA{$RD_valid}{$key}";
			}
		}
	}
	close INTER;close REF;
}elsif ($UserDB =~ /gene.sgRNA.db.fasta.txt$/) {
	print "# Reading User's fasta sequence database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/$label\_result_BO.txt");
	open (REF,">$dir/$label\_result_RO.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_sgRNA{$RD_valid}}){
			$seq_num{$seq}=0;
		}
	}
	while (<UDB>) {
		chomp;
		my ($inter_seq) = (split /\t/,$_)[-1];
		if (exists $seq_num{$inter_seq}) {
			$seq_num{$inter_seq}++;
		}
	}
	close UDB;
	print "# Find out the 'RO' and 'BO' sgRNA.\n";
	foreach my $RD_valid (@RD_valid) {
		foreach my $key (keys %{$RD_sgRNA{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				print INTER "$RD_sgRNA{$RD_valid}{$key}";
			}else{
				print REF "$RD_sgRNA{$RD_valid}{$key}";
			}
		}
	}
	close INTER;close REF;
}
print "# Your job is done, check result.\n# END\n\n";
my $oo = time() - $time;
print "Total time consumption is $oo second.\nDone!\n";
sub err {
	my ( $msg ) = @_;
	print "\n$msg\n\n";
	die "\n";
}

