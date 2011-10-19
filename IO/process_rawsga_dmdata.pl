use strict;

my %query_hash;

my $ind=1;
if(scalar(@ARGV) < 1) { print "usage: perl process_rawsga_dmdata.pl <filename>\n"; exit; }
my $currfile = $ARGV[0];


open TMP, "cut -f1 $currfile |";
while(<TMP>) {
    my $q = $_;
    chomp $q;

    if(not exists $query_hash{$q}) {
	$query_hash{$q}=$ind++;
    }
}
close TMP;

open TMP, "cut -f2 $currfile |";
while(<TMP>) {
    my $q = $_;
    chomp $q;

    if(not exists $query_hash{$q}) {
	$query_hash{$q}=$ind++;
    }
}
close TMP;


open IN,"< $currfile";
open OUT,"> $currfile"."_numeric";

while(<IN>) {
    my $line = $_;
    chomp $line;
    my @fields =  split(/\t/,$line);

    #grab query,array ids, and relevant fields
    #print OUT "$query_hash{$fields[0]}\t$query_hash{$fields[1]}\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7]\n";
    print OUT "$query_hash{$fields[0]}\t$query_hash{$fields[1]}\t".join("\t",@fields[2..(scalar(@fields)-1)])."\n";
}

close IN;
close OUT;

open OUT,">$currfile"."_orfidmap";
my @queries = keys %query_hash;
foreach my $q (@queries) {
    print OUT $q."\t$query_hash{$q}\n";
}
close OUT;

