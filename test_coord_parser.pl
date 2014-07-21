#!/usr/bin/env perl -w
use strict; 
use Moose; 
use Getopt::Long;
use Linkage::multiMI;
use Linkage::coord_parser;

my $window_size = 10;
my $location_file = "locations.in";
my $sites_file = "sites.in";
my $outfile_location = "mir.out";
GetOptions(
    "window|w=i" => \$window_size,
    "location|l=s" => \$location_file,
    "sites|s=s" => \$sites_file,
    "outfile|o=s" => \$outfile_location,
) or die "Usage: -w window_size -i snp.in -o mir.out";
die "Usage: ./test_coord_parser.pl -w window_size -l locations.in -s sites.in -o mir.out\n" if @ARGV != 0;

my $parser = Linkage::coord_parser->new(
    location_file => $location_file,
    sites_file => $sites_file,
    locations => [],
    site_info => []
);

my $data = Linkage::multiMI->new(
    region_1 => {}, # hash containing haplotype => frequency for region 1
    region_2 => {}, # hash containing haplotype => frequency for region 1
    joint_haplotypes => {} # hash with haplotype => frequency for joint region
);

$parser -> parseLocations(); # read data into memory
$parser -> parseSites();
my $sample_size = ${$parser -> {"site_info"}}[0];
print "Sample size: $sample_size \n";
my $tot_region_length = ${$parser -> {"location_info"}}[1];
print "Total region length: $tot_region_length\n";
open(TABFILE, ">$outfile_location") || die "Can't open output file $!";
print TABFILE "start_pos_1\tstart_pos_2\twindow_size\tmir\tsys_error\n";
for (my $i = 1; $i <= $tot_region_length + 1 - $window_size; $i+=$window_size) {
    $data -> {"region_1"} = {$parser -> getRegion($i, $window_size)};
    for (my $j = $i; $j < $tot_region_length + 1 - $window_size; $j+=$window_size) {
        $data -> {"region_2"} = {$parser -> getRegion($j, $window_size)};
        $data -> {"joint_haplotypes"} = 
            {$parser -> setJointRegion($i, $j, $window_size)};
        print TABFILE "$i\t$j\t$window_size\t", $data->mir(), "\t", $data->systematic_error($sample_size), "\n";
    }
}
close TABFILE;

exit;


