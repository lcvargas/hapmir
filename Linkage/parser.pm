package Linkage::parser;
use Moose; 

has 'snp_location' => ( isa => 'Any', is => 'rw',  required => 1);

sub openRegion1 {
    my $self = shift;
    my $region1_location = shift;
    my $window_size = shift;
    open(REGION1, $region1_location) || die "can't open file: $!";
    my %r1_uniques;
    while (<REGION1>) {
        chomp;
        my $line = $_;
        my $found_hap = 0; 
        foreach my $haplotype (keys %r1_uniques) {
            if ($haplotype =~ m/$line/) {
                $r1_uniques{$haplotype}++;
                $found_hap = 1;
                last;
            }
        }
        if ($found_hap == 0) {
            $r1_uniques{$line} = 1;
        }
    }
    close REGION1;
    print "region 1:\n";
    foreach my $haplotype (keys %r1_uniques) {
        print "$haplotype: $r1_uniques{$haplotype}\n";
    }
    return %r1_uniques;
}

sub openRegion2 {
    my $self = shift;
    my $region2_location = shift;
    open(REGION2, $region2_location) || die "can't open file: $!";
    my %r2_uniques;
    while (<REGION2>) {
        chomp;
        my $line = $_;
        my $found_hap = 0; 
        foreach my $haplotype (keys %r2_uniques) {
            if ($haplotype =~ m/$line/) {
                $r2_uniques{$haplotype}++;
                $found_hap = 1;
            last;
        }
        }
        if ($found_hap == 0) {
            $r2_uniques{$line} = 1;
        }
    }
    close REGION2;
    print "region 2:\n";
    foreach my $haplotype (keys %r2_uniques) {
        print "$haplotype: $r2_uniques{$haplotype}\n";
    }
    return %r2_uniques;
}

sub openJointRegion {
    my $self = shift;
    my $joint_region_location = shift;
    open(JOINTREGION, $joint_region_location) || die "can't open file: $!";
    my %joint_region_uniques;
    while (<JOINTREGION>) {
        chomp;
        my $line = $_;
        my $found_hap = 0; 
        foreach my $haplotype (keys %joint_region_uniques) {
            if ($haplotype =~ m/$line/) {
                $joint_region_uniques{$haplotype}++;
                $found_hap = 1;
                last;
            }
        }
        if ($found_hap == 0) {
            $joint_region_uniques{$line} = 1;
        }
    }
    close JOINTREGION;
    print "joint region:\n";
    foreach my $haplotype (keys %joint_region_uniques) {
        print "$haplotype: $joint_region_uniques{$haplotype}\n";
    }
    return %joint_region_uniques;
}

1;
