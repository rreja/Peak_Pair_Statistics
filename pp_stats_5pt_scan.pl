use strict;
use warnings;
use List::Util qw(max);
use Getopt::Std;


my %opt;
getopt('disgp',\%opt);
&help_message if (defined $opt{h});


my (%idx_files,%hash_O);
my (%pp_files,%hash_idx);
my $gsize;
my $pt_cutoff;
if (defined $opt{p}){
    $pt_cutoff = $opt{'p'};
}
else{
    $pt_cutoff = 5;
}

my $genome = $opt{'g'};
if($genome eq "NA"){
    $gsize = $opt{'s'};
}
else
{
   
    $gsize = find_gsize($genome);   
}


# opening the directory where all the index files are present.
my $idxdir = $opt{'i'};
$idxdir = check_dir($idxdir);
opendir IDX, $idxdir || die "Cannot open the directory";
while( (my $filename = readdir(IDX))){
    if(($filename =~ /.*idx.*/) || ($filename =~ /.*tab.*/)){
        my @fpart = split(/\./,$filename);
        my $key_name = shift(@fpart);
        $idx_files{$key_name} = 1;
        open IN,$idxdir.$filename || die "File not found";
        my $sum = 0;
        while(<IN>){
            chomp($_);
            next if($_ =~ /^#/);
            next if($_ =~ /^chrom/);
            my @cols = split(/\t/,$_);
            $sum = $sum + $cols[2] + $cols[3];
        }
        $hash_idx{$key_name} = $sum;
        close(IN);
     }
}
closedir(IDX);

# opening the directory where all the S_ and D_ files exist
my $dir = $opt{'d'};
$dir = check_dir($dir);
opendir DIR, $dir || die "Cannot open the directory";
while(my $filename = readdir(DIR)){
    if(($filename =~ /^S_/) || ($filename =~ /^O_/)){
        while(my($k,$v) = each %idx_files){
            if($filename =~ /$k/){ # Important step..Matching the idx file to the corresponding peak and orphan files.
                $pp_files{$filename} = $k;
            }
            
        }
        
    }
}
closedir(DIR);

open OUT, ">".$dir."peak_pair_stats.txt" || die "File not found";
print OUT "Filename\tPeak-pair_mode\tPeaks_in_peak-pairs\tOrphan_peaks\tMedian_peak-pair_occupancy\tAverage_peak-pair_occupancy\tFraction_of_all_mapped_reads_in_peak-pairs\n";
open OUT1, ">".$dir."signal2noise_qt_scan.txt" || die "File not found";
#print OUT "top_5pt_signal:noise\ttop_10pt_signal:noise\ttop_25pt_signal:noise\ttop_50pt_signal:noise\ttop_75pt_signal:noise\ttop_100pt_signal:noise\n";


while(my ($name, $idx) = each %pp_files)
{
    if($name =~ /^O_/){
        open IN, $dir.$name || die "File not found";
        my $count_O = 0;
        while(<IN>){
            chomp($_);
            next if(($_ =~ /^chrom/) || ($_ =~ /^#/));
            next if($_ =~ /^#/);
            $count_O++;
        }
        $hash_O{$idx} = $count_O;
        #print $count_O."\n";
        close(IN);
        
    }
}
while(my($name,$idx) = each %pp_files){
    if($name =~ /^S_/){
        open IN, $dir.$name || die "File not found";
        my @tmp = split(/_/,$name);
        my $count_S = 0;
        my $sum_of_col6 = 0;
        my (@cwdist,@tags,$D);
        if(pop(@tmp) =~ /^\D\d+\D(\d+).*/) # S_* files should have the same format sigma-exclusion-filter==> s5e20F1
        {$D = $1;}
        while(<IN>){
            chomp($_);
            next if($_ =~ /^#/);
            my @cols = split(/\t/,$_);
            push(@cwdist,(split(/=/,$cols[8]))[1]);
            push(@tags,$cols[5]);
            $count_S++;
            $sum_of_col6+=$cols[5];
                       
        }
        my $no_of_peaks = 2*$count_S;
        my $noise_per_bp = ($hash_idx{$idx} - $sum_of_col6)/($gsize - ($D*$count_S));
        my @sorted_tags = sort { $b <=> $a } @tags;
        my %quant = sum_quantile(\@sorted_tags,$pt_cutoff,$D,$noise_per_bp);
            
        my $vec = $idx."\t".mode(@cwdist)."\t".$no_of_peaks."\t".$hash_O{$pp_files{$name}}."\t".median(@tags)."\t".mean(@tags);
        #my $quants = $quant_1."\t".$quant_5."\t".$quant_10."\t".$quant_25."\t".$quant_50."\t".$quant_75."\t".$quant_100;
        print OUT $vec."\t".$sum_of_col6/$hash_idx{$idx}."\n";
        print OUT1 $idx."\n";
        print OUT1 "Quantile_range\tsignal2noise\n";
        for my $k ( sort {$a<=>$b} keys %quant){
            print OUT1 $k."\t".$quant{$k}."\n";
	    
        }
	
        close(IN);
    }
    
}

sub sum_quantile{
    my ($array,$cf,$ex,$noise) = @_;
    my $sum = 0;
    my $length = scalar(@$array);
    my %results = ();
    my $old_quant = 0;
    my $diff;
    for(my $cutoff=$cf;$cutoff<=100;$cutoff+=$cf){
        my $new_quant = int($length*($cutoff/100)); # this depdends on cutoff. It needs to be changed.
        $diff = $new_quant - $old_quant;
        for(;$old_quant<$new_quant;$old_quant++){
            $sum += $$array[$old_quant];
        }
    
    my $sum_per_bp = $sum/($ex*($diff));
    $results{$cutoff} = $sum_per_bp/$noise;
    #print $sum_per_bp/$noise ."\n";
    $sum = 0;
        
    }
    return %results;
}

sub median{
    my @a = sort {$a <=> $b} @_;
  return ($a[$#a/2] + $a[@a/2]) / 2;
}

sub mean{
  @_ or return 0;
  my $sum = 0;
  $sum += $_ foreach @_;
  return $sum/@_;
}

sub mode
{
    my %c;
    foreach my $e (@_) {
	$c{$e}++;
    }
    my $best = max(values %c);
    foreach my $n (keys %c){
        if($best == $c{$n}){
            return $n;
        }
    }
    
}

sub find_gsize
{
    my $genome = shift;
    if($genome eq "sg07"){
        
        return 15000000;
    }
    elsif($genome eq "mm09"){
        return(3000000000);
    }
    elsif($genome eq "mm08"){
        return(3000000000);
    }
    elsif($genome eq "hg18"){
        return(3000000000);
    }
    elsif($genome eq "hg19"){
        return(3000000000);
    }
    elsif($genome eq "dm03"){
        return(139500000);
    }
    
    
}

sub check_dir{
    my $dir = shift;
    my $tmp = substr($dir,-1,1);
    if($tmp eq "/"){
        return($dir);
    }
    else{
        return($dir."/");
    }
}

sub help_message {
  print qq{
Program: pp_stats_5pt_scan.pl (Calculate stats on peak pairs)
Contact: Rohit Reja <rzr142\@psu.edu>
Usage:   pp_stats_5pt_scan.pl -i <index_file_directory> -d <path_to_directory_containing_S_*_and_D_*_files> -g <organism_name> -p 5

    NOTE:    If you input files were saved using MS excel, then use:  perl -p -e 's/^M//g;' <input_file> > <input_file_no_excel_characters>
             to remove the excel characters in your file. ^M sould be typed as "ctrl-v-m". Or else the script will not work properly.
             
    Options: -i <path1>     path to the folder with index files[accepted index file extensions, idx, tab]. 
             -d <path2>     path to the folder with S_*.gff and O_* files. 
             -g             organism, sg07=>yeast, mm09=>MouseV9, mm08=>MouseV8, hg18=>human18, hg19=>human19, dm03=>Drosophila
             -s             size of genome[optional] In case of other genomes, set -g as NA and -s as the size of genome (see ex. below)
	     -p <number>    the percent quantile you need to use to scan. For ex. scan window of 5 is default.

    Example:
      perl pp_stats_5pt_scan.pl -i  ./ -d genetrack_s5e10F1/cwpair_output_mode_f0u0d100b3/ -g sg07
      perl pp_stats_5pt_scan.pl -i  ./ -d genetrack_s5e10F1/cwpair_output_mode_f0u0d100b3/ -g NA -s 160000000
      perl pp_stats_5pt_scan.pl -i  ./ -d genetrack_s5e10F1/cwpair_output_mode_f0u0d100b3/ -g NA -s 160000000 -p 10
      
    Output:
    Produces a "peak_pair_stats.txt" file in the folder that contains the S_*and D_* files. It also produces a
    "signal2noise_qt_scan.txt" file in the folder that contains the S_* files.
    
  
  };
  exit;
}


