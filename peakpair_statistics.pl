use strict;
use warnings;
use List::Util qw(max);

# opening the directory where all the S_ and D_ files exists
my $dir = $ARGV[0];  # remember to put a "/" at the end of the directory.
opendir my($dh), $dir || die "Cannot open the directory";
my @files = readdir $dh;
closedir $dh;

my @test = (2,3,2,2,2,2,4,5,5,5,5);
my %file;

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
        $file{$fname_parts[1]} = $count_O;
        #print $count_O."\n";
        close(IN);
        
    }
    elsif($fname_parts[0] eq "S"){
        open IN, $dir.$name || die "File not found";
        my $count_S = 0;
        my $sum_of_col6 = 0;
        my (@cwdist,@tags);
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
        my $vec = mode(@cwdist)."\t".$no_of_peaks."\t".median(@tags)."\t".mean(@tags)."\t".$sum_of_col6;
        print $vec."\n";
        close(IN);
    }
    
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

