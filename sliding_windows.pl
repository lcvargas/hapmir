#!/opt/local/bin/perl -w 
use strict;
use Moose;
use Linkage::multiMI;
use Linkage::parser;

my $parser = Linkage::parser->new(
    snp_location => "joint_region.snp"
);

my $data = Linkage::multiMI->new(
    region_1 => {},
    region_2 => {},
    joint_haplotypes => {}
);

die "Usage: ./sliding_windows.pl window_size outfile.csv\n" if @ARGV < 2;
my $window_size = $ARGV[0];
my $tot_region_length = $parser -> snpLength();
my %total_mir;
my %count_mir;
print "|| Window size || Reg1 start || Reg2 start ||     MIR     ||     Error     || \n";
for (my $i = 0; $i < $tot_region_length + 1 - $window_size; $i++) {
    for (my $j = 0; $j < $tot_region_length + 1 - $window_size; $j++) {
        $data -> {"region_1"} = {$parser -> setRegion1($i, $window_size)};
        $data -> {"region_2"} = {$parser -> setRegion2($j, $window_size)};
        $data -> {"joint_haplotypes"} = 
            {$parser -> setJointRegion($i, $j, $window_size)};
        print "||     $window_size      ||     $i      ||     $j      ||  " , 
            $data -> mir(), "     ||    ", 
            $data -> systematic_error(1000), "     || \n";
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

open(RFILE, "+>>$ARGV[1]") || die "can't open file $!";
seek(RFILE, 0, 0);
if (<RFILE> !~ m/^\"Window/) {
    print RFILE '"Window_size","Distance", "MIR"', "\n"; 
}
seek(RFILE, 0, 2);
for my $distance (keys %total_mir) {
    print RFILE "$window_size,$distance,", $total_mir{$distance} / $count_mir{$distance}, " \n";
}
close RFILE;


exit;


