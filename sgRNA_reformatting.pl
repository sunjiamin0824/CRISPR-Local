#!/usr/bin/perl
use strict;
use warnings;
die  "Usage:perl $0 <Designing_mode> <GFF_index> <All_match_file> <Outfile>\n" unless (@ARGV == 4);

open (INDEX,"<$ARGV[1]") or die;
my %exon_list_gene;
while (<INDEX>) {
	chomp;
	my($gene,$exon)=split /\t/,$_;
	$exon_list_gene{$exon}=$gene;
}
close INDEX;
open (SEQMAP,"<$ARGV[2]") or die;
readline SEQMAP;
open (OUTPUT,">$ARGV[3]") or die;
my $sgRNA_score;my $sgRNA_pam;my $sgRNA_pos;my $sgRNA_info;
my($OT_gene,$OT_pos,$OT_pam);
my %sgRNA_POS;my %sgRNA;
while (<SEQMAP>) {
	chomp;
	my($OT_id,undef,$OT_seq,$sgRNA_id,$sgRNA_seq,$mismatch) = split("\t",$_); 
	my $Exon = $sgRNA_id;
	$Exon =~ s/\#\[.+//;
	if (exists $exon_list_gene{$Exon}) {
		($OT_gene,$OT_pos,$OT_pam) = split ("\#",$OT_id);
		($sgRNA_info,$sgRNA_pos,$sgRNA_pam,$sgRNA_score) = (split "\#",$sgRNA_id)[1,2,3,4];
		if ($ARGV[0] eq "cpf1") {
			$OT_seq = $OT_pam.$OT_seq;
			$sgRNA_seq = $sgRNA_pam.$sgRNA_seq;
		}else{
			$OT_seq = $OT_seq.$OT_pam;
			$sgRNA_seq = $sgRNA_seq.$sgRNA_pam;
		}
		$sgRNA_POS{"$exon_list_gene{$Exon}:$sgRNA_pos"}{$Exon}=$sgRNA_info;
		$sgRNA{"$exon_list_gene{$Exon}:$sgRNA_pos"}{$OT_pos}="$exon_list_gene{$Exon}\t$sgRNA_pos\t$sgRNA_seq\t$sgRNA_score\t$OT_gene\t$OT_pos\t$OT_seq\t$mismatch";
	}
}
close SEQMAP;
my %Exon_list;my %Position_list;
foreach my $key1 (sort keys %sgRNA_POS) {
	foreach	my $key2 (sort keys %{$sgRNA_POS{$key1}}){
		$Exon_list{$key1}.="$key2;";
		$Position_list{$key1}.="$sgRNA_POS{$key1}{$key2};";
	}
	foreach my $key3 (sort keys %{$sgRNA{$key1}}) {
		print OUTPUT "$sgRNA{$key1}{$key3}\t$Exon_list{$key1}\t$Position_list{$key1}\n";
	}
}
close OUTPUT;
system("rm $ARGV[2]");
