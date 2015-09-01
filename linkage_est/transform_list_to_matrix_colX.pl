#!/usr/bin/perl

use strict;

if ($#ARGV != 1 ) {
	print "usage: transform_list_to_matrix.pl <path_to_list.txt> <column_to_use_PERL_MODE>\\n";
	exit;
}
my $path=$ARGV[0];
my @path_parts = split(/\//, $path);
my @files;
my $file = $path_parts[@path_parts-1];
my $folder;
if($file =~ /.txt/) { 
	pop(@path_parts); 
	$folder = join("/", @path_parts) . "/";
	$files[0] = $folder . $file;
} else { 
	$folder = join("/", @path_parts);
	@files = <$folder/*>;
}
my $col = $ARGV[1];
print "Folder:\t" . $folder ."\n";
print "Files:\t" . join("\n\t", @files) . "\n";
print "Column (PERL mode):\t" . $col . "\n";

# LOAD the ORF coordinates
open(COORD, "/Users/Anastasia/Laboratory/Datasets/Utils_SGA/chrom_coordinates.txt");
my %orf2start;
while(<COORD>) {
	my $line = $_;
	chomp($line);
	my @fields = split(/\t/, $line);
	$orf2start{$fields[0]} = &min($fields[2], $fields[3]);
}
close COORD;


for(my $f = 0; $f < @files; $f++) {

	my $j = $f + 1;
	print "File $j ... ";

	open(BP_FILE, $files[$f]);
	open(OUTPUT, '>' . $files[$f] . '_matrix_col' . $col . '.txt');
	
	my $count_q=0;
	my $count_a=0;
	my @matrix;
	
	my %query_to_id;
	my %array_to_id;
	
	# CHECK!!! VERY IMPORTANT!!!
	my $make_symmetric = 0;

	while(<BP_FILE>) {    
		my $line = $_;
		chomp($line);
		my @fields = split(/\t/, $line);
		
		#unless ($fields[0] =~ /^Y/ && $fields[1] =~ /^Y/) { next; }
	 
		my $query = $fields[0];
		my $array = $fields[1];
		#print $query . "\t" . $array . "\n";
		  
		if ($query eq 'undefined' | $array eq 'undefined' | $array eq 'YOR202W') { next; }
		
		unless(exists $query_to_id{$query}) {
			$query_to_id{$query} = $count_q;
			$count_q++;
		}
		
		unless(exists $array_to_id{$array}) {
			$array_to_id{$array} = $count_a;
			$count_a++;
		}
		   
		$matrix[$query_to_id{$query}][$array_to_id{$array}] = $fields[$col];
		
		# Comment out when the matrix should not be symmetric
		if($make_symmetric == 1) {
			unless(exists $query_to_id{$array}) {
				$query_to_id{$array} = $count_q;
				$count_q++;
			}
		
			unless(exists $array_to_id{$query}) {
				$array_to_id{$query} = $count_a;
				$count_a++;
			}
			
			$matrix[$query_to_id{$array}][$array_to_id{$query}] = $fields[$col];		}
		#-----------------------
	
	}
	close BP_FILE;
	
	# SORT QUERIES AND ARRAYS
	my @queries = sort sort_by_chrom keys %query_to_id;
	#my @queries = keys %query_to_id;
	my @arrays = sort sort_by_chrom keys %array_to_id;
	
	
	print OUTPUT "ORF\tGWEIGHT\t";
	
	for(my $i = 0; $i < @arrays; $i++) {
		print OUTPUT $arrays[$i];

		if($i < @arrays - 1) {
			print OUTPUT "\t";
		} else {
			print OUTPUT "\n";
		}
	}
	
	print OUTPUT "EWEIGHT\t\t";

	for(my $i = 0; $i < @arrays; $i++) {
		print OUTPUT 1;

		if($i < @arrays - 1) {
			print OUTPUT "\t";

		} else {
			print OUTPUT "\n";
		}
	}
	
	for(my $i = 0; $i < @queries; $i++) {
		print OUTPUT $queries[$i] . "\t1\t";

		my $query_ind = $query_to_id{$queries[$i]};
		for(my $j = 0; $j < @arrays; $j++) {
			my $array_ind = $array_to_id{$arrays[$j]};
			
			if( (length($matrix[$query_ind][$array_ind]) == 0) ) { print OUTPUT '0'; }
			else { print OUTPUT $matrix[$query_ind][$array_ind]; }
			
			if($j < @arrays - 1) {
				print OUTPUT "\t";
			} else {
				print OUTPUT "\n";
			}
		}
	}
	
	close OUTPUT;
	
	print "done\n";
}


sub min { 
	if ($_[0]>$_[1]) {return $_[1]} else {return $_[0]}; 
} 

sub sort_by_chrom {
	my @fields_a = split(/_/, $a);
	my @fields_b = split(/_/, $b);
	if(substr($a,1,1) eq substr($b,1,1)) { return $orf2start{$fields_a[0]} <=> $orf2start{$fields_b[0]}; }
	else { return substr($a,1,1) cmp substr($b,1,1); }
}
