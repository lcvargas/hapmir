package Linkage::parser;
use Moose; 

has 'snp_location' => ( isa => 'Any', is => 'rw',  required => 1);

sub snpLength {
    my $self = shift;
    open(DATA, $self->{"snp_location"}) || die "can't open file: $!";
    my $first_line = <DATA>;
    chomp($first_line);
    close DATA;
    return length($first_line);
}

# setRegion1($starting_position, $window_size)
sub setRegion1 {
    my $self = shift;
    my $starting_position = shift;
    my $window_size = shift;
    open(DATA, $self->{"snp_location"}) || die "can't open file: $!";
    my %r1_uniques;
    my $line;
    while (<DATA>) {
        chomp;
        $line = substr($_, $starting_position, $window_size);
        if (exists $r1_uniques{$line}) {
            $r1_uniques{$line}++;
        }
        else {
            $r1_uniques{$line} = 1;
        }
    }
    close DATA;
    # print "region 1:\n";
    # foreach my $haplotype (keys %r1_uniques) {
    #     print "$haplotype: $r1_uniques{$haplotype}\n";
    # }
    return %r1_uniques;
}

# setRegion2($starting_position, $window_size)
sub setRegion2 {
    my $self = shift;
    my $starting_position = shift;
    my $window_size = shift;
    open(DATA, $self->{"snp_location"}) || die "can't open file: $!";
    my %r2_uniques;
    my $line;
    while (<DATA>) {
        chomp;
        $line = substr($_, $starting_position, $window_size);
        if (exists $r2_uniques{$line}) {
            $r2_uniques{$line}++;
        }
        else {
            $r2_uniques{$line} = 1;
        }
    }
    close DATA;
    # print "region 2:\n";
    # foreach my $haplotype (keys %r2_uniques) {
    #     print "$haplotype: $r2_uniques{$haplotype}\n";
    #  }
    return %r2_uniques;
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
