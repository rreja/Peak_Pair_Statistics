use strict;
use warnings;

# opening the directory where all the S_ and D_ files exists
my $dir = $ARGV[0];  # remember to put a "/" at the end of the directory.
opendir my($dh), $dir || die "Cannot open the directory";
my @files = readdir $dh;
closedir $dh;

my %replicates;

foreach my $name (@files)
{
    next if(($name =~ /^\./) || ($name eq "statistics.txt"));
    my @fname_parts = split(/\_/,$name);
    if($fname_parts[0] eq "O"){
        open IN, $dir.$name || die "File not found";
        my $count_O = 0;
        while(<IN>){
            next if(($_ =~ /^chrom/) || ($_ =~ /^#/));
            $count_O++;
            
            
        }
        print $count_O."\n";
        
    }
    elsif($fname_parts[0] eq "S"){
        
    }
    
}


