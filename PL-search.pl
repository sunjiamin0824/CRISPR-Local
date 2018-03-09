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
$label ||= "PL-search";
my $errmsg = "Use '$0 -h' for quick help; for more information, please see README.";
my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,
	perl $0 -g ZmB73_paralogous_gene.list -i ZmB73.reference.database.txt -o /your_dir/

	or

	perl $0 -g ZmB73_paralogous_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/

	
--- Options ---

	-h: show this help and exit

	-g <string>	:Paralog gene list
	-i <string>	:Reference sgRNA database (RD)
	-u <string>	:User's sgRNA database (UD)
	-o <string>	:Output path (defualt: $dir_default)
	-l <label>	:Name prefix for output file (default:PL_search)
};
err( $helpmsg ) if $opt_help;

err( "ERROR - No paralog gene list file provided.\n$errmsg" ) if !$List;
err( "ERROR - No reference database file provided.\n$errmsg" ) if !$RefDB;

print "\n  Welcome to CRISPR-Local\n";
print " Program PL_search: a local tool for search exclusive and commom target for paralogous gene pair.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 1.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjm0824\@webmail.hzau.edu.cn"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today : $local_time\n\n";
print "# Program PL-search: a local tool for search exclusive and commom target for paralogous gene pair.\n";

print "# Reading paralog gene list file.\n";
open (LIST,$List) or die "Can't open $List for reading!\n";

my %gene_list;my %gene_pair;my(%RD_sgRNA,%RD_pos);

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
	my($gene,$pos,$seq,$score,$OT_1,$OT_2,$OT_3,$OT_4,$OT_score)=(split /\t/,$_)[0,1,2,3,4,5,6,7,-1];
	if (exists $gene_list{$gene}) {
		$RD_pos{"$gene:$seq"}.="$gene:$pos:$score;";
		$RD_sgRNA{$gene}{$seq}="$gene\t$pos\t$seq\t$score\t$OT_1\t$OT_2\t$OT_3\t$OT_4\t$OT_score";
	}
}
close SGRNA;

my %Inter_hash;
foreach my $gene_pair (sort keys %gene_pair) {
	my @pair = (split ",",$gene_pair);
	my @inter_seq = keys %{$RD_sgRNA{$pair[0]}};
	for (my $i=0;$i<=$#pair;$i++) {
		@inter_seq = grep {$RD_sgRNA{$pair[$i]}{$_}} @inter_seq;
	}
	foreach my $inter (@inter_seq) {
		$Inter_hash{$gene_pair}{$inter} = 1;
		foreach my $paralog (@pair) {
			delete $RD_sgRNA{$paralog}{$inter};
		}
	}
}
if (!$UserDB) {
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/$label\_common_targets.txt") or die "Can't open $label\_common_targets.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		my @inter_seq = keys %{$Inter_hash{$gene_pair}};
		foreach my $inter (@inter_seq) {
			print OUT "$inter";
			foreach my $paralog (@pair) {
				print OUT "\t".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;
	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/$label\_exclusive_targets.txt") or die "Can't open $label\_exclusive_targets.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_sgRNA{$pair[$i]}};
			for (my $j=0;$j<=$#RD_seq;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$RD_seq[$j]}\n";
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
	open (OUT,">$dir/$label\_common_targets.txt") or die "Can't open $label\_common_targets.txt for writing!\n";
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
				print OUT "\t".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/$label\_exclusive_targets.txt") or die "Can't open $label\_exclusive_targets.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_sgRNA{$pair[$i]}};
			my %UD;
			foreach my $RD_seq (@RD_seq) {
				if (exists $UD_score{$pair[$i]}{$RD_seq}) {
					$UD{$pair[$i]}{$RD_seq} = "UD\t$UD_score{$pair[$i]}{$RD_seq}";
				}else{
					$UD{$pair[$i]}{$RD_seq} = "RD\t0";
				}
			}
			for (my $j=0;$j<=$#RD_seq;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$RD_seq[$j]}\t$UD{$pair[$i]}{$RD_seq[$j]}\n";
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
	open (OUT,">$dir/$label\_common_targets.txt") or die "Can't open $label\_common_targets.txt for writing!\n";
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
				print OUT "\t".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/$label\_exclusive_targets.txt") or die "Can't open $label\_exclusive_targets.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_sgRNA{$pair[$i]}};
			my %UD;
			foreach my $RD_seq (@RD_seq) {
				if (exists $UD_score{$pair[$i]}{$RD_seq}) {
					$UD{$pair[$i]}{$RD_seq} = "UD\t$UD_score{$pair[$i]}{$RD_seq}";
				}else{
					$UD{$pair[$i]}{$RD_seq} = "RD\t0";
				}
			}
			for (my $j=0;$j<=$#RD_seq;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$RD_seq[$j]}\t$UD{$pair[$i]}{$RD_seq[$j]}\n";
			}
		}
	}
	close OUT;
}elsif ($UserDB =~ /gene.sgRNA.db.fasta.txt/) {
	print "# Reading User's database file.\n";
	open (UDB,$UserDB) or die;
	while (<UDB>) {
		chomp;
		my($num,$seq) = (split /\t/,$_)[-3,-1];
		$UD_score{$seq}++;
	}
	close UDB;
	print "# Finding the commom targets...\n";
	open (OUT,">$dir/$label\_common_targets.txt") or die "Can't open $label\_common_targets.txt for writing!\n";
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
				print OUT "\t".$RD_pos{"$paralog:$inter"};
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;

	print "# Finding the exclusive targets...\n";
	open (OUT,">$dir/$label\_exclusive_targets.txt") or die "Can't open $label\_exclusive_targets.txt for writing!\n";
	foreach my $gene_pair (sort keys %gene_pair) {
		my @pair = (split ",",$gene_pair);
		for (my $i=0;$i<=$#pair;$i++) {
			my @RD_seq = sort keys %{$RD_sgRNA{$pair[$i]}};
			my %UD;
			foreach my $RD_seq (@RD_seq) {
				if (exists $UD_score{$pair[$i]}{$RD_seq}) {
					$UD{$pair[$i]}{$RD_seq} = "UD\t$UD_score{$pair[$i]}{$RD_seq}";
				}else{
					$UD{$pair[$i]}{$RD_seq} = "RD\t0";
				}
			}
			for (my $j=0;$j<=$#RD_seq;$j++) {
				print OUT "$RD_sgRNA{$pair[$i]}{$RD_seq[$j]}\t$UD{$pair[$i]}{$RD_seq[$j]}\n";
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
