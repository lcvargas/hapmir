#!/usr/bin/env perl -w
use strict; use Moose; use Getopt::Long;
use Linkage::multiMI;
use Linkage::parser;

my $window_size = 10;
my $snp_location = "snp.in";
my $outfile_location = "mir.out";
GetOptions(
    "window|w=i" => \$window_size,
    "snp|i|infile=s" => \$snp_location,
    "outfile|o=s" => \$outfile_location
) or die "Usage: -w window_size -i snp.in -o mir.out";
die "Usage: ./sliding_windows.pl -w window_size -i snp.in -o mir.out\n" if @ARGV != 0;

my $parser = Linkage::parser->new(
    snp_location => $snp_location,
    data => {}
);

my $data = Linkage::multiMI->new(
    region_1 => {}, # hash containing haplotype => frequency for region 1
    region_2 => {}, # hash containing haplotype => frequency for region 1
    joint_haplotypes => {} # hash with haplotype => frequency for joint region
);

my $tot_region_length = $parser -> snpLength();
my $sample_size = $parser -> parse(); # read data into memory and 
                                    # return the number of lines in file
# print "$sample_size\n";
open(TABFILE, ">$outfile_location") || die "Can't open output file $!";
print TABFILE "start_pos_1\tstart_pos_2\twindow_size\tmir\tsys_error\n";
for (my $i = 0; $i < $tot_region_length + 1 - $window_size; $i+=$window_size) {
    $data -> {"region_1"} = {$parser -> getRegion($i, $window_size)};
    for (my $j = $i; $j < $tot_region_length + 1 - $window_size; $j+=$window_size) {
        $data -> {"region_2"} = {$parser -> getRegion($j, $window_size)};
        $data -> {"joint_haplotypes"} = 
            {$parser -> setJointRegion($i, $j, $window_size)};
        print TABFILE "$i\t$j\t$window_size\t", $data->mir(), "\t", $data->systematic_error($sample_size), "\n";
        # produces mir_data.tab 
        # if (exists $total_mir{abs($i - $j)}) {
        #    $total_mir{abs($i - $j)} += $data -> mir();
        #    $total_error{abs($i - $j)} += $data ->     
        #        systematic_error($sample_size);
        #    $count{abs($i - $j)}++; 
            
        # }
        # else {
        #    $total_mir{abs($i - $j)} = $data -> mir();
        #     $total_error{abs($i - $j)} = $data -> 
        #        systematic_error($sample_size);
        #     $count{abs($i - $j)} = 1; 
        # }
    }
    # print "Completed: $i / ", $tot_region_length + 1 - $window_size, "\n";
}
close TABFILE;


# open(RFILE, "+>>$outfile_location") || die "can't open output file $!";
# seek(RFILE, 0, 0);
# if (<RFILE> !~ m/^Window/) {
#     print RFILE "Window_size\tDistance\tMIR\tError\n"; 
# }
# seek(RFILE, 0, 2);
# for my $distance (keys %total_mir) {
#    print RFILE "$window_size\t$distance\t", $total_mir{$distance} /  $count{$distance}, "\t", $total_error{$distance} / $count{$distance}, "\n";
# }
# close RFILE;


exit;


