#!/opt/local/bin/perl -w 
use strict;
use Moose;
use Linkage::multiMI;
use Linkage::parser;

my $parser = Linkage::parser->new(
    snp_location => "joint_region.snp"
);

my $data = Linkage::multiMI->new(
    region_1 => {$parser -> setRegion1(0, 10)},
    region_2 => {$parser -> setRegion2(10, 10)},
    joint_haplotypes => {$parser -> setJointRegion(0, 10, 10)}
);

print "snp length: " , $parser -> snpLength(), "\n";
print "Entropy at 1: ", $data -> entropy("region_1"), "\n"; # entropy at locus 1
print "Entropy at 2: ", $data -> entropy("region_2"), "\n"; # entropy at locus 2
print "Mutual information: ", $data -> mutual_information(), "\n";
print "MIR: ", $data -> mir(), "\n";
print "Systematic error: ", $data -> systematic_error(1000), "\n";
print "T statistic: ", $data -> t_statistic(1000), "\n";
exit;


