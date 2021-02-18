#!/usr/bin/perl

use strict;
use warnings;


open (IN,"$ARGV[0]") || die;   ###chr1.gene/chr1.promoter

my %hash;
while(<IN>){
   chomp;
   my @aa=split(/\s+/,$_);
     $hash{$aa[-1]}[0]++;     ###total C in this region;
     $hash{$aa[-1]}[1]++ if($aa[5] eq "CG");  ###total CG in this region;
     $hash{$aa[-1]}[2]++ if($aa[5] eq "CHG");  ###total CHG in this region;
     $hash{$aa[-1]}[3]++ if($aa[5] eq "CHH");  ###total CHH in this region;

     $hash{$aa[-1]}[4]++ if($aa[3]+$aa[4] >=4);  ###total covered C in this region;
     $hash{$aa[-1]}[5]++ if($aa[3]+$aa[4] >=4 && $aa[5] eq "CG");  ###total covered CG in this region;
     $hash{$aa[-1]}[6]++ if($aa[3]+$aa[4] >=4 && $aa[5] eq "CHG");  ###total covered CHG in this region;   
     $hash{$aa[-1]}[7]++ if($aa[3]+$aa[4] >=4 && $aa[5] eq "CHH");  ###total covered CHH in this region;

     $hash{$aa[-1]}[8] +=$aa[3] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]);   ###mReads num in C;
     $hash{$aa[-1]}[9] +=$aa[4] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]);  ###umReads num in C;

     $hash{$aa[-1]}[10] +=$aa[3] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CG");# mReads num in CG
     $hash{$aa[-1]}[11] +=$aa[4] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CG");# umReads num in CG : song_yunjie:20180322

     $hash{$aa[-1]}[12] +=$aa[3] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CHG");# mReads num in CHG
     $hash{$aa[-1]}[13] +=$aa[4] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CHG");# umReads num in CHG

     $hash{$aa[-1]}[14] +=$aa[3] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CHH");# mReads num in CHH
     $hash{$aa[-1]}[15] +=$aa[4] if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7] && $aa[5] eq "CHH");# umReads num in CHH

     $hash{$aa[-1]}[16] ++ if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]); ###mC num in this region;
     $hash{$aa[-1]}[17] ++ if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CG"); ###mCG num in this region;
     $hash{$aa[-1]}[18] ++ if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CHG"); ###mCHG num in this region;
     $hash{$aa[-1]}[19] ++ if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CHH"); ###mCHH num in this region;
     
     $hash{$aa[-1]}[20] +=$aa[3]/($aa[3]+$aa[4]) if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]); ###mC sum methylation level;
     $hash{$aa[-1]}[21] +=$aa[3]/($aa[3]+$aa[4]) if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CG"); ###mCG sum methylation level;
     $hash{$aa[-1]}[22] +=$aa[3]/($aa[3]+$aa[4]) if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CHG"); ###mCHG sum methylation level;
     $hash{$aa[-1]}[23] +=$aa[3]/($aa[3]+$aa[4]) if($aa[3]+$aa[4] >=4 && $aa[3] >=$aa[7]&& $aa[5] eq "CHH"); ###mCHH sum methylation level;
     
     $hash{$aa[-1]}[24] +=$aa[3]+$aa[4] if($aa[3]+$aa[4] >=4);    ###total reads  of covered C);
     $hash{$aa[-1]}[25] +=$aa[3]+$aa[4] if($aa[3]+$aa[4] >=4&& $aa[5] eq "CG");   ###total reads  of covered CG);
     $hash{$aa[-1]}[26] +=$aa[3]+$aa[4] if($aa[3]+$aa[4] >=4&& $aa[5] eq "CHG");   ###total reads  of covered CHG);
     $hash{$aa[-1]}[27] +=$aa[3]+$aa[4] if($aa[3]+$aa[4] >=4&& $aa[5] eq "CHH");   ###total reads  of covered CHG);
}

close(IN);



#print "Gene_id\tTotal_C\tCovered_C_num\tmReads_C_num\tumReads_C_num\tC_Methy_level\tTotal_CG\tCovered_CG_num\tmReads_CG_num\tumReads_CG_num\tCG_Methy_level\tTotal_CHH\tCovered_CHH_num\tmReads_CHH_num\tumReads_CHH_num\tCHH_Methy_level\tTotal_CHG\tCovered_CHG_num\tmReads_CHG_num\tumReads_CHG_num\tCHG_Methy_level\tTotal_methylated_C_num\n";

print "Gene_id\tTotal_C\tTotal_CG\tTotal_CHG\tTotal_CHH\tTotal_covered_C\tTotal_covered_CG\tTotal_covered_CHG\tTotal_covered_CHH\tmReads_C_num\tumReads_C_num\tmReads_CG_num\tumReads_CG_num\tmReads_CHG_num\tumReads_CHG_num\tmReads_CHH_num\tumReads_CHH_num\tTotal_methy_C_num\tTotal_methy_CG_num\tTotal_methy_CHG_num\tTotal_methy_CHH_num\tC_Methy_level_sum\tCG_Methy_level_sum\tCHG_Methy_level_sum\tCHH_Methy_level_sum\ttotal_reads_of_covered_C\ttotal_reads_of_covered_CG\ttotal_reads_of_covered_CHG\ttotal_reads_of_covered_CHH\tC_Methy_level\tCG_Methy_level\tCHG_Methy_level\tCHH_Methy_level\tC_Methy_density\tCG_Methy_density\tCHG_Methy_density\tCHH_Methy_density\tC_Methy_level_average\tCG_Methy_level_average\tCHG_Methy_level_average\tCHH_Methy_level_average\n";

#foreach my $kk(sort {(split(/\:/,$a))[1] <=> (split(/\:/,$b))[1]} keys %hash){
foreach my $kk( keys %hash){
    for(my $i=0;$i<=27;$i++){if(!$hash{$kk}[$i]){$hash{$kk}[$i]=0;}}
	my $m_level=0;
        my $m_level_CG=0;
	my $m_level_CHH=0;
	my $m_level_CHG=0;
	my $m_density=0;
	my $m_density_CG=0;
	my $m_density_CHG=0;
	my $m_density_CHH=0;
	my $m_average=0;
	my $m_average_CG=0;
	my $m_average_CHG=0;
	my $m_average_CHH=0;
    if($hash{$kk}[24]>0){
	   $m_level=sprintf("%.2f",100*$hash{$kk}[8]/($hash{$kk}[24]));
	}
    if($hash{$kk}[25]>0){
           $m_level_CG=sprintf("%.2f",100*$hash{$kk}[10]/($hash{$kk}[25]));
        }
    if($hash{$kk}[26]>0){
           $m_level_CHG=sprintf("%.2f",100*$hash{$kk}[12]/($hash{$kk}[26]));
        }
    if($hash{$kk}[27]>0){
           $m_level_CHH=sprintf("%.2f",100*$hash{$kk}[14]/($hash{$kk}[27]));
        }

    if($hash{$kk}[4]>0){
           $m_density=sprintf("%.2f",100*$hash{$kk}[16]/($hash{$kk}[4]));
        }
    if($hash{$kk}[5]>0){
           $m_density_CG=sprintf("%.2f",100*$hash{$kk}[17]/($hash{$kk}[5]));
        }
    if($hash{$kk}[6]>0){
           $m_density_CHG=sprintf("%.2f",100*$hash{$kk}[18]/($hash{$kk}[6]));
        }    
    if($hash{$kk}[7]>0){
           $m_density_CHH=sprintf("%.2f",100*$hash{$kk}[19]/($hash{$kk}[7]));
        }
     
    if($hash{$kk}[4]>0){
           $m_average=sprintf("%.2f",100*$hash{$kk}[20]/($hash{$kk}[4]));
        }

    if($hash{$kk}[5]>0){
           $m_average_CG=sprintf("%.2f",100*$hash{$kk}[21]/($hash{$kk}[5]));
        }
    if($hash{$kk}[6]>0){
           $m_average_CHG=sprintf("%.2f",100*$hash{$kk}[22]/($hash{$kk}[6]));
        }
    if($hash{$kk}[7]>0){
           $m_average_CHH=sprintf("%.2f",100*$hash{$kk}[23]/($hash{$kk}[7]));
        }
	$hash{$kk}[28]=$m_level;
        $hash{$kk}[29]=$m_level_CG;
	$hash{$kk}[30]=$m_level_CHG;
        $hash{$kk}[31]=$m_level_CHH;
	$hash{$kk}[32]=$m_density;
	$hash{$kk}[33]=$m_density_CG;
	$hash{$kk}[34]=$m_density_CHG;
	$hash{$kk}[35]=$m_density_CHH;
	$hash{$kk}[36]=$m_average;
	$hash{$kk}[37]=$m_average_CG;
	$hash{$kk}[38]=$m_average_CHG;
	$hash{$kk}[39]=$m_average_CHH;
	#print $kk;
	my @tem=@{$hash{$kk}};
	#print $kk;
	print $kk,"\t";
	#print join("\t",(split(/\:/,$kk))),"\t";
        #print join("\t",($tem[0],$tem[2],$tem[4],$tem[5],$tem[17],$tem[1],$tem[3],$tem[7],$tem[8],$tem[18],$tem[9],$tem[10],$tem[11],$tem[12],$tem[19],$tem[13],$tem[14],$tem[15],$tem[16],$tem[20],$tem[6])),"\n";
	print join("\t",@tem),"\n";
}


