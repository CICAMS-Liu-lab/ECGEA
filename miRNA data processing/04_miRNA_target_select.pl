#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION="1.0";
GetOptions( \%opts,"t=s","n=i","i=s","head=s","h!");


my $usage = <<"USAGE";
        Usage:perl $0 -i index_list -t file -n index_col
	-i	indexfile
	-t	tablefile
	-n	index line num in tablefile
	-head	T or F
	-h      Display this usage information

USAGE

die $usage if ( !( $opts{i} && $opts{t} && $opts{n} ) || $opts{h} );


my $index_file = $opts{i};
my $table_file = $opts{t};
my $index_col = $opts{n}-1;
my $head = $opts{head}?$opts{head}:"T";
open (FIN1, "<$index_file");
my %index;
while(<FIN1>){
	chomp;
	$index{$_}=1;
	}
open (FIN2, "<$table_file");
if($head eq "T"){
	my $head_line = <FIN2>;
	print $head_line;
	}
while(<FIN2>){
         chomp;    
         my @line =split(/\t+/,$_); 
         if (exists $index{$line[$index_col]}){
		print $_."\n";
		}else{
	}
}
close FIN1;
close FIN2;
