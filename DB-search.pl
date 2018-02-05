#!/usr/bin/perl
use strict;
use warnings;

print "Hello, World...\n";
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

my ($List, $RefDB, $UserDB, $Path, $opt_help,$opt_number);
my $time = time();

GetOptions(	
			'l=s' => \$List,
			'i=s' => \$RefDB,
			'u:s' => \$UserDB,
			'o=s' => \$Path,
			'h!'  => \$opt_help,		#Help message
			'N:i' => \$opt_number,
          );

my $dir_default = getcwd;             #default output
$Path ||= $dir_default;
my $dir =$Path;
my $errmsg = "Use '$0 -h' for quick help; for more information, please see README.";

my $helpmsg = qq{
=============== CRISPR-Local ===============

--- Usage ---
e.g.,
	perl $0 -l ZmB73_query_gene.list -i ZmB73.reference.database.txt -o /your_dir/ -N 3
	
	or

	perl $0 -l ZmB73_query_gene.list -i ZmB73.reference.database.txt -u ZmC01.gene.sgRNA.db.alignment.txt -o /your_dir/ -N 3

--- Options ---

	-h: show this help and exit

	-l: <Query_gene_list>
	-i: <Reference_database>
	-u: <User's database>
	-o: <Output_path> (defualt: $dir_default)
	-N: <int> Report up to top <int> targets per gene (default: all result)

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

err( "ERROR - No query gene list file provided.\n$errmsg" ) if !$List;
err( "ERROR - No reference database file provided.\n$errmsg" ) if !$RefDB;


if ($opt_number) {
	err ("ERROR - Wrong option of number of targets to reproted : -N $opt_number (supported: > 0).\n$errmsg") if ($opt_number <= 0);
}

print "\n  Welcome to CRISPR-Local\n";
print " Program DB_search: a local tool for search best target of query gene.\n";
print "  ---------------------------------------------------------\n";
print " Version   : 1.0"."\n";
print " Copyright : Free software"."\n";
print " Author    : Jiamin Sun"."\n";
print " Email     : sunjiamin0824\@qq.com"."\n\n";

my $local_time;
$local_time = localtime();
print "# Today : $local_time\n\n";
print "# Program DB-search: Search the reference and user's database to select the sgRNA with high on-target score, low off-target score.\n";
my $ran = int(rand(100000));
print "\n# Your qry_id is $ran.\n";
print "# Reading gene list file.\n";
open (LIST,$List) or die "Can't open $List for reading!\n";
my %gene_list;my(%RD_sgRNA,%UD_sgRNA,%RD_score,%UD_score);my(%RD_score_num,%RD_sgRNA_num);my %hash;
while (<LIST>) {
	chomp;
	$gene_list{$_}=undef;
}
my @list = keys %gene_list;
my $list = @list;
print "# There are $list gene in your list.\n";
print "# Reading reference database file.\n";
open (SGRNA,$RefDB) or die $!;
while (<SGRNA>) {
	chomp;
	my($gene,$pos,$seq,$score,$OT_gene,$OT_pos,$OT_seq,$mismatch,$CFD)=(split /\t/,$_)[0,1,2,3,4,5,6,7,-1];
	if (exists $gene_list{$gene}) {
		if (exists $RD_sgRNA_num{$gene}{$score} and exists $RD_score_num{$gene}{$seq}) {
			$RD_sgRNA{$gene}{$score}.="\n$gene\t$pos\t$score\t$seq\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch\t$CFD";
		}elsif (exists $RD_sgRNA_num{$gene}{$score}) {
			$RD_sgRNA{$gene}{$score}.="\n$gene\t$pos\t$score\t$seq\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch\t$CFD";
			$RD_score{$gene}{$seq}=$score;
			$RD_score_num{$gene}{$seq} = 1;
		}elsif (exists $RD_score_num{$gene}{$seq}) {
			$RD_sgRNA{$gene}{$score}="$gene\t$pos\t$score\t$seq\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch\t$CFD";
			$RD_sgRNA_num{$gene}{$score} = 1;
			$RD_score{$gene}{$seq}.=",$score";
			$RD_score_num{$gene}{$seq}++;
		}else{
			$RD_sgRNA{$gene}{$score}="$gene\t$pos\t$score\t$seq\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch\t$CFD";
			$RD_sgRNA_num{$gene}{$score} = 1;
			$RD_score{$gene}{$seq}=$score;
			$RD_score_num{$gene}{$seq} = 1;
		}
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
	open (LEFT,">$dir/Invalid_gene_RD_$ran.list") or die "Can't open Invalid_gene_RD_$ran.list for writing!\n";
	foreach my $key (sort keys %gene_left) {
		print LEFT "$key\n";
	}
	print "# The invalid gene is output to Invalid_gene_RD_$ran.list, please add annotation information of these invalid gene.\n";
	close LEFT;
}
if (!$UserDB) {
	print "# Output reference database search result.\n";
	open (OUT,">$dir/Gene_search_result_RD_$ran.txt") or die "Can't open Gene_search_result_RD_$ran.txt for writing!\n";
	foreach my $RD_valid (@RD_valid) {
		my @rs2 = sort desc_sort_subject(keys %{$RD_sgRNA{$RD_valid}});
		if (!$opt_number) {
			for (my $i=0;$i<=$#rs2;$i++) {
				print OUT "$RD_sgRNA{$RD_valid}{$rs2[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($rs2[$i]) {
					print OUT "$RD_sgRNA{$RD_valid}{$rs2[$i]}\n";
				}
			}
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
			$UD_sgRNA{$gene}{$pos}="$_\n";
			$UD_score{$gene}{$seq}.="$pos,";
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
		open (LEFT,">$dir/Invalid_gene_UD_$ran.list") or die "Can't open Invalid_gene_UD_$ran.list for writing!\n";
		foreach my $key (sort keys %gene_left) {
			print LEFT "$key\n";
		}
		close LEFT;
		print "# The invalid gene is output to Invalid_gene_UD_$ran.list.\n";
	}
	print "# Find out the 'reference-database-only', 'user-database-only' and 'Common' sgRNA.\n";
	my @inter = grep {$RD_sgRNA{$_}} @UD_valid;
	my %merge = map {$_ => 1} @RD_valid,@UD_valid;
	my @merge = sort keys (%merge);
	my @RDC_gene = grep {!$RD_sgRNA{$_}} @merge;
	my @UDC_gene = grep {!$UD_sgRNA{$_}} @merge;
	open (REF,">$dir/Gene_search_result_RD_only_$ran.txt");
	open (USER,">$dir/Gene_search_result_UD_only_$ran.txt");
	open (CORNA,">$dir/Gene_search_result_Co-sgRNA_$ran.txt");
	print "# Output the 'RD', 'UD' and 'Common' sgRNA result.\n";
	foreach my $RDC_gene (@RDC_gene) {
		my @rs2 = sort keys %{$UD_sgRNA{$RDC_gene}};
		if (!$opt_number) {
			for (my $i=0;$i<=$#rs2;$i++) {
				print USER "$UD_sgRNA{$RDC_gene}{$rs2[$i]}";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($rs2[$i]) {
					print USER "$UD_sgRNA{$RDC_gene}{$rs2[$i]}";
				}
			}
		}
	}
	foreach my $UDC_gene (@UDC_gene) {
		my @rs2 = sort desc_sort_subject(keys %{$RD_sgRNA{$UDC_gene}});
		if (!$opt_number) {
			for (my $i=0;$i<=$#rs2;$i++) {
				print REF "$RD_sgRNA{$UDC_gene}{$rs2[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($rs2[$i]) {
					print REF "$RD_sgRNA{$UDC_gene}{$rs2[$i]}\n";
				}
			}
		}
	}
	foreach my $inter (sort @inter) {
		my @RD_seq = sort keys %{$RD_score{$inter}};
		my @UD_seq = sort keys %{$UD_score{$inter}};
		my @inter_seq = grep {$RD_score{$inter}{$_}} @UD_seq;
		my %merge_seq = map {$_ => 1} @RD_seq,@UD_seq;
		my @merge_seq = sort keys (%merge_seq);
		my @RDC_seq = grep {!$RD_score{$inter}{$_}} @merge_seq;
		my @UDC_seq = grep {!$UD_score{$inter}{$_}} @merge_seq;
		my(@score_RDC,@score_UDC,@score_inter);
		foreach my $RDC_seq (@RDC_seq) {
			my @row = split ",",$UD_score{$inter}{$RDC_seq};
			push(@score_RDC,@row);
		}
		foreach my $UDC_seq (@UDC_seq) {
			if ($RD_score_num{$inter}{$UDC_seq} > 1) {
				my @scores = split ",",$RD_score{$inter}{$UDC_seq};
				push(@score_UDC,@scores);
			}else{
				push(@score_UDC,$RD_score{$inter}{$UDC_seq});
			}
		}
		$hash{$_} = 1, foreach @score_UDC;
		@score_UDC = keys %hash;
		undef %hash;
		foreach my $inter_seq (@inter_seq) {
			if ($RD_score_num{$inter}{$inter_seq} > 1) {
				my @scores = split ",",$RD_score{$inter}{$inter_seq};
				push(@score_UDC,@scores);
			}else{
				push(@score_inter,$RD_score{$inter}{$inter_seq});
			}
		}
		$hash{$_} = 1, foreach @score_inter;
		@score_inter = keys %hash;
		undef %hash;
		@score_UDC = sort desc_sort_subject(@score_UDC);
		@score_inter = sort desc_sort_subject(@score_inter);
		if (!$opt_number) {
			for (my $i=0;$i<=$#score_RDC;$i++) {
				print USER "$UD_sgRNA{$inter}{$score_RDC[$i]}";
			}
			for (my $i=0;$i<=$#score_UDC;$i++) {
				print REF "$RD_sgRNA{$inter}{$score_UDC[$i]}\n";
			}
			for (my $i=0;$i<=$#score_inter;$i++) {
				print CORNA "$RD_sgRNA{$inter}{$score_inter[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($score_RDC[$i]) {
					print USER "$UD_sgRNA{$inter}{$score_RDC[$i]}";
				}
			}
			for (my $i=0;$i<$opt_number;$i++) {
				if ($score_UDC[$i]) {
					print REF "$RD_sgRNA{$inter}{$score_UDC[$i]}\n";
				}
			}
			for (my $i=0;$i<$opt_number;$i++) {
				if ($score_inter[$i]) {
					print CORNA "$RD_sgRNA{$inter}{$score_inter[$i]}\n";
				}
			}
		}
	}
	close REF;close USER;close CORNA;
}elsif ($UserDB =~ /intergenic.sgRNA.db.alignment.txt$/) {
	print "# Reading User's intergenic alignment database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/Gene_search_result_CO-sgRNA_$ran.txt");
	open (REF,">$dir/Gene_search_result_RD_only_$ran.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_score{$RD_valid}}){
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
	print "# Find out the 'reference-database-only' and 'Common' sgRNA.\n";
	my @seq_num = keys %seq_num;
	foreach my $RD_valid (@RD_valid) {
		my(@co_num,@RD_num);
		foreach my $key (keys %{$RD_score{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@co_num,@scores);
				}else{
					push(@co_num,$RD_score{$RD_valid}{$key});
					$RD_sgRNA{$RD_valid}{$RD_score{$RD_valid}{$key}}.="\t$seq_num{$key}";
				}
			}else{
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@RD_num,@scores);
				}else{
					push(@RD_num,$RD_score{$RD_valid}{$key});
				}
			}
		}
		$hash{$_} = 1, foreach @co_num;
		@co_num = keys %hash;
		undef %hash;

		$hash{$_} = 1, foreach @RD_num;
		@RD_num = keys %hash;
		undef %hash;

		@co_num = sort desc_sort_subject(@co_num);
		@RD_num = sort desc_sort_subject(@RD_num);
		if (!$opt_number) {
			for (my $i=0;$i<=$#co_num;$i++) {
				print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
			}
			for (my $i=0;$i<=$#RD_num;$i++) {
				print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($co_num[$i]) {
					print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
				}
			}
			for (my $i=0;$i<$opt_number;$i++) {
				if ($RD_num[$i]) {
					print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
				}
			}
		}
	}
	close INTER;close REF;
}elsif ($UserDB =~ /gene.sgRNA.db.fastq.txt$/) {
	print "# Reading User's fastq reads database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/Gene_search_result_CO-sgRNA_$ran.txt");
	open (REF,">$dir/Gene_search_result_RD_only_$ran.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_score{$RD_valid}}){
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
	print "# Find out the 'reference-database-only' and 'Common' sgRNA.\n";
	my @seq_num = keys %seq_num;
	foreach my $RD_valid (@RD_valid) {
		my(@co_num,@RD_num);
		foreach my $key (keys %{$RD_score{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@co_num,@scores);
				}else{
					push(@co_num,$RD_score{$RD_valid}{$key});
					$RD_sgRNA{$RD_valid}{$RD_score{$RD_valid}{$key}}.="\t$seq_num{$key}";
				}
			}else{
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@RD_num,@scores);
				}else{
					push(@RD_num,$RD_score{$RD_valid}{$key});
				}
			}
		}
		$hash{$_} = 1, foreach @co_num;
		@co_num = keys %hash;
		undef %hash;

		$hash{$_} = 1, foreach @RD_num;
		@RD_num = keys %hash;
		undef %hash;

		@co_num = sort desc_sort_subject(@co_num);
		@RD_num = sort desc_sort_subject(@RD_num);
		if (!$opt_number) {
			for (my $i=0;$i<=$#co_num;$i++) {
				print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
			}
			for (my $i=0;$i<=$#RD_num;$i++) {
				print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($co_num[$i]) {
					print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
				}
			}
			for (my $i=0;$i<$opt_number;$i++) {
				if ($RD_num[$i]) {
					print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
				}
			}
		}
	}
	close INTER;close REF;
}elsif ($UserDB =~ /gene.sgRNA.db.fasta.txt$/) {
	print "# Reading User's fasta sequence database file.\n";
	my %seq_num;
	open (UDB,$UserDB) or die;
	open (INTER,">$dir/Gene_search_result_CO-sgRNA_$ran.txt");
	open (REF,">$dir/Gene_search_result_RD_only_$ran.txt");
	foreach my $RD_valid (@RD_valid) {
		foreach my $seq (keys %{$RD_score{$RD_valid}}){
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
	print "# Find out the 'reference-database-only' and 'Common' sgRNA.\n";
	my @seq_num = keys %seq_num;
	foreach my $RD_valid (@RD_valid) {
		my(@co_num,@RD_num);
		foreach my $key (keys %{$RD_score{$RD_valid}}) {
			if ($seq_num{$key} > 0) {
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@co_num,@scores);
				}else{
					push(@co_num,$RD_score{$RD_valid}{$key});
					$RD_sgRNA{$RD_valid}{$RD_score{$RD_valid}{$key}}.="\t$seq_num{$key}";
				}
			}else{
				if ($RD_score_num{$RD_valid}{$key} > 1) {
					my @scores = split ",",$RD_score{$RD_valid}{$key};
					push(@RD_num,@scores);
				}else{
					push(@RD_num,$RD_score{$RD_valid}{$key});
				}
			}
		}
		$hash{$_} = 1, foreach @co_num;
		@co_num = keys %hash;
		undef %hash;

		$hash{$_} = 1, foreach @RD_num;
		@RD_num = keys %hash;
		undef %hash;

		@co_num = sort desc_sort_subject(@co_num);
		@RD_num = sort desc_sort_subject(@RD_num);
		if (!$opt_number) {
			for (my $i=0;$i<=$#co_num;$i++) {
				print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
			}
			for (my $i=0;$i<=$#RD_num;$i++) {
				print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
			}
		}else{
			for (my $i=0;$i<$opt_number;$i++) {
				if ($co_num[$i]) {
					print INTER "$RD_sgRNA{$RD_valid}{$co_num[$i]}\n";
				}
			}
			for (my $i=0;$i<$opt_number;$i++) {
				if ($RD_num[$i]) {
					print REF "$RD_sgRNA{$RD_valid}{$RD_num[$i]}\n";
				}
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
sub desc_sort_subject {
	$b <=> $a;
}