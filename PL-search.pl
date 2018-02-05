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

my ($List, $RefDB, $UserDB, $Path, $opt_help,);
my $time = time();
GetOptions(	
			'l=s' => \$List,
			'i=s' => \$RefDB,
			'u:s' => \$UserDB,
			'o=s' => \$Path,
			'h!'  => \$opt_help,		#Help message
          );

my $dir_default = getcwd;             #default output
$Path ||= $dir_default;
my $dir =$Path;
my $errmsg = "Use '$0 -h' for quick help; for more information, please see README.";

my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,
	perl $0 -l ZmB73_paralogous_gene.list -i ZmB73.reference.database.txt -o /your_dir/

	or

	perl $0 -l ZmB73_paralogous_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/ -N 3

	
--- Options ---

	-h: show this help and exit

	-l: <paralog_gene_list>
	-i: <Reference_database>
	-u: <User's database>
	-o: <Output_path> (defualt: $dir_default)

};
err( $helpmsg ) if $opt_help;

if (not ($List && $RefDB) ) {
	print
		"\n" .
		"Please input the option, and press <Enter>.\n" .
		"For quick help, input '-h' and press <Enter>.\n" .
		"\n";
	die "\n";
}

err( "ERROR - No paralog gene list file provided.\n$errmsg" ) if !$List;
err( "ERROR - No reference database file provided.\n$errmsg" ) if !$RefDB;


print "\n  Welcome to CRISPR-Local\n";
print " Program PL_search: a local tool for search exclusive and commom target for paralogous gene pair.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 1.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjiamin0824\@qq.com"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today : $local_time\n\n";
print "# Program PL-search: a local tool for search exclusive and commom target for paralogous gene pair.\n";
my $ran = int(rand(100000));
print "\n# Your mission id is $ran.\n";

print "# Reading paralog gene list file.\n";
open (LIST,$List) or die "Can't open $List for reading!\n";
my %gene_list;my %gene_pair;my(%RD_sgRNA,%RD_score,%RD_pos);

while (<LIST>) {
	chomp;
	my @gene = (split ",",$_);
	foreach my $gene (@gene) {
		$gene_list{$gene} = undef;
	}
	$gene_pair{$_} = undef;
}

print "# Reading reference database file.\n";
open (SGRNA,$RefDB) or die $!;
while (<SGRNA>) {
	chomp;
	my($gene,$pos,$seq,$score,$OT_gene,$OT_pos,$OT_seq,$mismatch,$CFD)=(split /\t/,$_)[0,1,2,3,4,5,6,7,-1];
	if (exists $gene_list{$gene}) {
		$RD_pos{"$gene:$seq"}="$pos:$score";
		$RD_sgRNA{$gene}{$score}="$gene\t$pos\t$score\t$seq\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch\t$CFD";
		$RD_score{$gene}{$seq}=$score;
	}
}
close SGRNA;

my %Inter_hash;
foreach my $gene_pair (sort keys %gene_pair) {
	my @pair = (split ",",$gene_pair);
	my @inter_seq = keys %{$RD_score{$pair[0]}};
	for (my $i=0;$i<=$#pair;$i++) {
		@inter_seq = grep {$RD_score{$pair[$i]}{$_}} @inter_seq;
	}
	foreach my $inter (@inter_seq) {
		$Inter_hash{$gene_pair}{$inter} = 1;
		foreach my $paralog (@pair) {
			delete $RD_score{$paralog}{$inter};
		}
	}
}

if (!$UserDB) {
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/Paralog_search_result_common_$ran.txt") or die "Can't open Paralog_search_result_common_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		my @inter_seq = keys %{$Inter_hash{$gene_pair}};
		foreach my $inter (@inter_seq) {
			print OUT "$inter";
			foreach my $paralog (@pair) {
				print OUT "\t$paralog:".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;
	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/Paralog_search_result_exclusive_$ran.txt") or die "Can't open Paralog_search_result_exclusive_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_score{$pair[$i]}};
			my @score_RD;
			foreach my $RD_seq (@RD_seq) {
				push(@score_RD,$RD_score{$pair[$i]}{$RD_seq});
			}
			@score_RD = sort desc_sort_subject(@score_RD);
			for (my $j=0;$j<=$#score_RD;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$score_RD[$j]}\n";
			}
		}
	}
	close OUT;
	print "# Your job is done, check result.\n# END\n\n";
	my $oo = time() - $time;
	print "Total time consumption is $oo second.\nDone!\n";
	exit(0);
}

my %UD_merge;my %UD_score;
if ($UserDB =~ /gene.sgRNA.db.alignment.txt$/) {
	print "# Reading User's database file.\n";
	open (UDB,$UserDB) or die;
	while (<UDB>) {
		chomp;
		my($gene,$num,$seq) = (split /\t/,$_)[0,2,-1];
		if (exists $gene_list{$gene}) {
			$UD_score{$gene}{$seq}+=$num;
		}
	}
	close UDB;
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/Paralog_search_result_common_$ran.txt") or die "Can't open Paralog_search_result_common_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		my @inter_seq = keys %{$Inter_hash{$gene_pair}};
		for (my $i=0;$i<=$#pair;$i++) {
			my @UD_seq = keys %{$UD_score{$pair[$i]}};
			$UD_merge{$_} = 1, foreach @UD_seq;
		}
		foreach my $inter (@inter_seq) {
			print OUT "$inter";
			if (exists $UD_merge{$inter}) {
				print OUT "\tUD";
			}else{
				print OUT "\tRD";
			}
			foreach my $paralog (@pair) {
				print OUT "\t$paralog:".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/Paralog_search_result_exclusive_$ran.txt") or die "Can't open Paralog_search_result_exclusive_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_score{$pair[$i]}};
			my @score_RD;my %UD;
			foreach my $RD_seq (@RD_seq) {
				push(@score_RD,$RD_score{$pair[$i]}{$RD_seq});
				if (exists $UD_score{$pair[$i]}{$RD_seq}) {
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "UD\t$UD_score{$pair[$i]}{$RD_seq}";
				}else{
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "RD\t0";
				}
			}
			@score_RD = sort desc_sort_subject(@score_RD);
			for (my $j=0;$j<=$#score_RD;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$score_RD[$j]}\t$UD{$score_RD[$j]}\n";
			}
		}
	}
	close OUT;
}elsif ($UserDB =~ /(intergenic.sgRNA.db.alignment.txt)|(gene.sgRNA.db.fastq.txt)/) {
	print "# Reading User's database file.\n";
	open (UDB,$UserDB) or die;
	while (<UDB>) {
		chomp;
		my($num,$seq) = (split /\t/,$_)[-3,-1];
		$UD_score{$seq}+=$num;
	}
	close UDB;
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/Paralog_search_result_common_$ran.txt") or die "Can't open Paralog_search_result_common_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		my @inter_seq = keys %{$Inter_hash{$gene_pair}};
		foreach my $inter (@inter_seq) {
			print OUT "$inter";
			if (exists $UD_score{$inter}) {
				print OUT "\tUD";
			}else{
				print OUT "\tRD";
			}
			foreach my $paralog (@pair) {
				print OUT "\t$paralog:".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/Paralog_search_result_exclusive_$ran.txt") or die "Can't open Paralog_search_result_exclusive_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_score{$pair[$i]}};
			my @score_RD;my %UD;
			foreach my $RD_seq (@RD_seq) {
				push(@score_RD,$RD_score{$pair[$i]}{$RD_seq});
				if (exists $UD_score{$RD_seq}) {
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "UD\t$UD_score{$RD_seq}";
				}else{
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "RD\t0";
				}
			}
			@score_RD = sort desc_sort_subject(@score_RD);
			for (my $j=0;$j<=$#score_RD;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$score_RD[$j]}\t$UD{$score_RD[$j]}\n";
			}
		}
	}
	close OUT;
}elsif ($UserDB =~ /gene.sgRNA.db.fasta.txt/) {
	print "# Reading User's database file.\n";
	open (UDB,$UserDB) or die;
	while (<UDB>) {
		chomp;
		my($seq) = (split /\t/,$_)[-1];
		$UD_score{$seq}++;
	}
	close UDB;
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/Paralog_search_result_common_$ran.txt") or die "Can't open Paralog_search_result_common_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		my @inter_seq = keys %{$Inter_hash{$gene_pair}};
		foreach my $inter (@inter_seq) {
			print OUT "$inter";
			if (exists $UD_score{$inter}) {
				print OUT "\tUD";
			}else{
				print OUT "\tRD";
			}
			foreach my $paralog (@pair) {
				print OUT "\t$paralog:".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/Paralog_search_result_exclusive_$ran.txt") or die "Can't open Paralog_search_result_exclusive_$ran.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_score{$pair[$i]}};
			my @score_RD;my %UD;
			foreach my $RD_seq (@RD_seq) {
				push(@score_RD,$RD_score{$pair[$i]}{$RD_seq});
				if (exists $UD_score{$RD_seq}) {
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "UD\t$UD_score{$RD_seq}";
				}else{
					$UD{$RD_score{$pair[$i]}{$RD_seq}} = "RD\t0";
				}
			}
			@score_RD = sort desc_sort_subject(@score_RD);
			for (my $j=0;$j<=$#score_RD;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$score_RD[$j]}\t$UD{$score_RD[$j]}\n";
			}
		}
	}
	close OUT;
}

print "# Your job is done, check result.\n# END\n\n";
my $oo = time() - $time;
print "Total time consumption is $oo second.\nDone!\n";

sub err {
	my ( $msg ) = @_;
	print "\n$msg\n\n";
	die "\n";
}
sub desc_sort_subject {
	$b <=> $a;
}
