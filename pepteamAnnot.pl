#!/usr/bin/perl

use strict;
use warnings;
use Text::CSV;
use Getopt::Long;
use File::Basename;
    

my $csv= Text::CSV->new();
##ajouter le get_opt..
my $inputCSV="";#=$ARGV[0];##the file containing all proteome id from ensembl+severals attributes related to this protein
my $inputTSV="";#=$ARGV[1];##the Scoring output TSV file (tab separated value) containing three columns Ensembl protein id/P-value/Z-score

my $species;#=$ARGV[3];##the species 
my @row;
my %rows_table;
my @header_table;
my $prog = basename $0;
my $format="";
my $sep;
my $outputCSV="";


GetOptions( 'input1=s' => \$inputCSV, 'input2=s' => \$inputTSV, 'output=s' => \$outputCSV);
if ( $inputTSV eq "" or $inputCSV eq "") {
    print "usage: ",$prog," --input1=TARGETFILE --input2=QUERYFILE --output=OUTPUTFILE [-verbose] [--help]\n";
    exit;
}
if($format eq ""){
	$format="csv";
}
if($outputCSV eq ""){
	$outputCSV=	"./out1.$format";
}
if($format eq "csv"){
	$sep=",";
}
else{
	$sep="\t";

}
#############################################################
#load the protein database and index it
#############################################################
print STDERR "loading protein annotations...\n";
my $columns_number;
my $counter=0;

open (CSV, "<", $inputCSV) or die $!;
while (<CSV>) {
	#start parsing fields
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$columns_number=$#columns;
	
		if($columns[0] ne "" && $columns[0] ne "predicted gene"){## ensembl protein id not null and not a predicted gene (to check if it's necessary)
			
			for(my $i=1;$i<=$columns_number;$i++){
				push(@row,$columns[$i]);

			}
			$rows_table{$columns[0]}= [@row];
	
			$counter++;
			##display some progress bar
			if($counter%10000==0){
				#print STDERR ".";
			}
			@row=();
		}
			
	}
	else {
		my $err = $csv->error_input;
		print STDERR "Failed to parse line: $err";
	}
}
close CSV;

#############################################################
#remove the header count
#############################################################
$counter=$counter-1;
#############################################################
##display the total number of valids proteins loaded
#############################################################

print STDERR "   ...$counter proteins loaded\n";

#############################################################
# open query TSV file and retrieve related informations 
#############################################################
open (TSV, "<", $inputTSV) or die $!;
#############################################################
# open target TSV and CSV file  
#############################################################
open(OUTPUT,">$outputCSV") or die ("Error in creating newFile") ;
#open(CSVOUTPUT,">$outputCSV") or die ("Error in creating newFile") ;

my $id_query;
my $id_target;
my @split;

#############################################################
# for each Id in Query TSV file, use hashtable to retrieve 
# related informations
#############################################################
print STDERR "annotating proteins...\n";
$counter=0;
while (my $ligne = <TSV>) {

	if ($ligne=~/(.*)\s(.*)\s(.*)\n/){
		
		foreach my $Id  (keys(%rows_table)){
			if($Id eq $1){
				@row=$rows_table{$Id};
				print OUTPUT "\"$Id\"$sep\"$2\"$sep\"$3\"$sep";
				#print CSVOUTPUT "\"$Id\",\"$2\",\"$3\",";
				
				for(my $i=0;$i<$columns_number;$i++){
			
					if($i!=$columns_number-1){
						#print TSVOUTPUT "\"$rows_table{$Id}[$i]\"\t";
						print OUTPUT "\"$rows_table{$Id}[$i]\"$sep";
					}
					else{
						print OUTPUT "\"$rows_table{$Id}[$i]\"\n";
						#print CSVOUTPUT "\"$rows_table{$Id}[$i]\"\n";
					}
				}
				$counter++;			
			}
		}
		if($counter%100==0){
			#print STDERR ".";
		}		
	}
}
print STDERR "   ...annotated $counter proteins\n";
close(TSV);

close(OUTPUT);
#close(CSVOUTPUT);








