#!/opt/local/bin/perl -w 
use strict;
use Moose;
use Getopt::Long;
use Linkage::multiMI;
use Linkage::parser;

my $window_size;
my $snp_location;
my $outfile_location = "mir.out";
die "Usage: ./sliding_windows.pl -w window_size -i infile.snp -o outfile.csv\n" if @ARGV != 6;
GetOptions(
    "window|w=i" => \$window_size,
    "snp|i|infile=s" => \$snp_location,
    "outfile|o=s" => \$outfile_location
) or die "Usage: -w window_size -i infile.snp -o outfile.csv";

my $parser = Linkage::parser->new(
    snp_location => $snp_location 
);

my $data = Linkage::multiMI->new(
    region_1 => {},
    region_2 => {},
    joint_haplotypes => {}
);

my $tot_region_length = $parser -> snpLength();
my %total_mir;
my %count_mir;
for (my $i = 0; $i < $tot_region_length + 1 - $window_size; $i++) {
    for (my $j = 0; $j < $tot_region_length + 1 - $window_size; $j++) {
        $data -> {"region_1"} = {$parser -> setRegion1($i, $window_size)};
        $data -> {"region_2"} = {$parser -> setRegion2($j, $window_size)};
        $data -> {"joint_haplotypes"} = 
            {$parser -> setJointRegion($i, $j, $window_size)};
        if (exists $total_mir{abs($i - $j)}) {
            $total_mir{abs($i - $j)} += $data -> mir();
            $count_mir{abs($i - $j)}++; 
        }
        else {
            $total_mir{abs($i - $j)} = $data -> mir();
            $count_mir{abs($i - $j)} = 1; 
        }
    }
}

open(RFILE, "+>>$outfile_location") || die "can't open output file $!";
seek(RFILE, 0, 0);
if (<RFILE> !~ m/^Window/) {
    print RFILE "Window_size\tDistance\tMIR\n"; 
}
seek(RFILE, 0, 2);
for my $distance (keys %total_mir) {
    print RFILE "$window_size\t$distance\t", $total_mir{$distance} / $count_mir{$distance}, " \n";
}
close RFILE;


exit;


