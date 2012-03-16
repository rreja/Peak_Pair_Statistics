use strict;
use warnings;
use List::Util qw(max);
use Getopt::Std;

my %opt;
getopt('disg',\%opt);

# opening the directory where all the S_ and D_ files exists
my $dir = $opt{'d'};  # remember to put a "/" at the end of the directory.
opendir my($dh), $dir || die "Cannot open the directory";
my @files = readdir $dh;
closedir $dh;

# opening the directory where all the index files are present.
my $idxdir = $opt{'i'}; #remember to put a "/" at the end of the directory.
opendir my($idh), $idxdir || die "Cannot open the directory";
my @indexes = readdir $idh;
closedir $idh;

my $gsize;
my (%hash_idx, %hash_O,%hash_S);
my $genome = $opt{'g'};
if($genome eq "NA"){
    $gsize = $opt{'s'};
}
else
{
   
    $gsize = find_gsize($genome);   
}



foreach my $idx (@indexes){
    next if($idx =~ /^\./);
    next if($idx !~ /.*\..*/);
    my @name = split(/\./,$idx);
    if(($name[1] eq "idx") || ($name[1] eq "tab")){
        open IN,$idxdir.$idx || die "File not found";
        my $sum = 0;
        while(<IN>){
            chomp($_);
            next if($_ =~ /^chrom/);
            my @cols = split(/\t/,$_);
            $sum = $sum + $cols[2] + $cols[3];
        }
        $hash_idx{$name[0]} = $sum;
        close(IN);
        
    }
    
    
}

print "Filename\tPeak-pair-mode\tPeaks-in-peak-pairs\tOrphan-peaks\tMedian-peak-pair-occupancy\tAverage-peak-pair-occupancy\tFraction-of-all-mapped-reads-in-peak-pairs\ttop-1pt-signal:noise\t";
print "top-5pt-signal:noise\ttop-10pt-signal:noise\ttop-25pt-signal:noise\ttop-50pt-signal:noise\ttop-75pt-signal:noise\ttop-100pt-signal:noise\n";

foreach my $name (@files)
{
    next if(($name =~ /^\./) || ($name eq "statistics.txt"));
    my @fname_parts = split(/\_/,$name);
    if($fname_parts[0] eq "O"){
        open IN, $dir.$name || die "File not found";
        my $count_O = 0;
        while(<IN>){
            chomp($_);
            next if(($_ =~ /^chrom/) || ($_ =~ /^#/));
            $count_O++;
        }
        $hash_O{$fname_parts[1]} = $count_O;
        #print $count_O."\n";
        close(IN);
        
    }
    elsif($fname_parts[0] eq "S"){
        open IN, $dir.$name || die "File not found";
        my $count_S = 0;
        my $sum_of_col6 = 0;
        my (@cwdist,@tags,$D);
        if($fname_parts[2] =~ /^\D\d+\D(\d+).*/)
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
        my $noise_per_bp = ($hash_idx{$fname_parts[1]} - $sum_of_col6)/($gsize - ($D*$count_S));
        my @sorted_tags = sort { $b <=> $a } @tags;
        my $quant_1 = sum_quantile(\@sorted_tags,1,$D,$noise_per_bp);
        my $quant_5 = sum_quantile(\@sorted_tags,5,$D,$noise_per_bp);
        my $quant_10 = sum_quantile(\@sorted_tags,10,$D,$noise_per_bp);
        my $quant_25 = sum_quantile(\@sorted_tags,25,$D,$noise_per_bp);
        my $quant_50 = sum_quantile(\@sorted_tags,50,$D,$noise_per_bp);
        my $quant_75 = sum_quantile(\@sorted_tags,75,$D,$noise_per_bp);
        my $quant_100 = sum_quantile(\@sorted_tags,100,$D,$noise_per_bp);
        
        my $vec = $fname_parts[1]."\t".mode(@cwdist)."\t".$no_of_peaks."\t".$hash_O{$fname_parts[1]}."\t".median(@tags)."\t".mean(@tags);
        my $quants = $quant_1."\t".$quant_5."\t".$quant_10."\t".$quant_25."\t".$quant_50."\t".$quant_75."\t".$quant_100;
        print $vec."\t".$sum_of_col6/$hash_idx{$fname_parts[1]}."\t".$quants."\n";
        close(IN);
    }
    
}

sub sum_quantile{
    my ($array,$cutoff,$ex,$noise) = @_;
    my $sum = 0;
    my $length = scalar(@$array);
    my $quant = int($length*($cutoff/100));
    for(my $i=0;$i<$quant;$i++){
        $sum += $$array[$i];
    }
    my $sum_per_bp = $sum/($ex*$quant);
    return($sum_per_bp/$noise);
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
    if($genome eq "sg7"){
        
        return 15000000;
    }
    elsif($genome eq "mm9"){
        return(3000000000);
    }
    elsif($genome eq "mm8"){
        return(3000000000);
    }
    elsif($genome eq "hg18"){
        return(3000000000);
    }
    elsif($genome eq "hg19"){
        return(3000000000);
    }
    elsif($genome eq "dm3"){
        return(120000000);
    }
    
    
}




