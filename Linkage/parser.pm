package Linkage::parser;
use Moose; 

has 'snp_location' => ( isa => 'Any', is => 'rw',  required => 1);
has 'data' => ( isa => 'ArrayRef', is => 'rw', required => 0);

sub snpLength {
    my $self = shift;
    open(DATA, $self->{"snp_location"}) || die "Can't open snp file: $!";
    my $first_line = <DATA>;
    chomp($first_line);
    close DATA;
    return length($first_line);
}

sub sampleSize {
    my $self = shift;
    open(DATA, $self->{"snp_location"}) || die "Can't open snp file: $!";
    while (<DATA>) {};
    return $.;

}

# getRegion($starting_position, $window_size)
sub getRegion {
    my $self = shift;
    my $starting_position = shift;
    my $window_size = shift;
    open(DATA, $self->{"snp_location"}) || die "can't open file: $!";
    my %uniques;
    my $line;
    while (<DATA>) {
        chomp;
        $line = substr($_, $starting_position, $window_size);
        if ($line !~ m/(0|1)+/) { # only allow 0's and 1's 
            next;
        }
        # throw out ?'d code here at this point; line must have 0 or 1 only; next if no match
        if (exists $uniques{$line}) {
            $uniques{$line}++;
        }
        else {
            $uniques{$line} = 1;
        }
    }
    close DATA;
    # print "region 1:\n";
    # foreach my $haplotype (keys %r1_uniques) {
    #     print "$haplotype: $r1_uniques{$haplotype}\n";
    # }
    return %uniques;
}


# setJointRegion($region_1_starting_position, $region_2_starting_position,
#   $window_size)
sub setJointRegion {
    my $self = shift;
    my $region_1_starting_position = shift;
    my $region_2_starting_position = shift;
    my $window_size = shift;
    open(DATA, $self->{"snp_location"}) || die "can't open file: $!";
    my %unique_haplotypes;
    my $line;
    while (<DATA>) {
        chomp;
        $line = substr($_, $region_1_starting_position, $window_size) . 
            substr($_, $region_2_starting_position, $window_size);
        if ($line !~ m/(0|1)+/) { # only allow 0's and 1's 
            next;
        }
        # throw out ?'d code here at this point; line must have 0 or 1 only; next if no match
        if (exists $unique_haplotypes{$line}) {
            $unique_haplotypes{$line}++;
        }
        else {
            $unique_haplotypes{$line} = 1;
        }
    }
    close DATA;
    # print "joint region:\n";
    # foreach my $haplotype (keys %unique_haplotypes) {
    #    print "$haplotype: $unique_haplotypes{$haplotype}\n";
    # }
    return %unique_haplotypes;
}


1;
