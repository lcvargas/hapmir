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

die "Usage: ./sliding_windows.pl window_size \n" if @ARGV < 1;
my $window_size = $ARGV[0];
my $tot_region_length = $parser -> snpLength();
my $region_1;
my $region_2;
my %MIRs;
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
    }
}
exit;


