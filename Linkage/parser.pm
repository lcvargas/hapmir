package Linkage::parser;
use Moose; 

has 'snp_location' => ( isa => 'Any', is => 'rw',  required => 1);
has 'data' => ( isa => 'HashRef', is => 'rw',  required => 0);

sub snpLength {
    my $self = shift;
    open(DATA, $self->{"snp_location"}) || die "Can't open snp file: $!";
    my $first_line = <DATA>;
    chomp($first_line);
    close DATA;
    return length($first_line);
}

# parse() reads the data into memory and returns the number of lines in the file
sub parse {
    my $self = shift;
    my $counter = 0;
    open(DATA, $self->{"snp_location"}) || die "Can't open snp file: $!";
    while (<DATA>) {
        chomp;
        $counter++;
        ${$self->{"data"}}{$counter} = $_;
    }
    close DATA;
    # foreach my $keys (keys %{$self->{"data"}}) {
    #    print "$keys => ${$self->{'data'}}{$keys}\n";
    # }
    return $counter; 
}

# getRegion($starting_position, $window_size)
sub getRegion {
    my $self = shift;
    my $starting_position = shift;
    my $window_size = shift;
    my %uniques;
    my $region;
    foreach my $line (keys %{$self->{"data"}}) {
        $region = substr(${$self->{'data'}}{$line}, $starting_position, $window_size);
        if ($region !~ m/(0|1)+/) { # only allow 0's and 1's 
            next;
        }
        if (exists $uniques{$region}) {
            $uniques{$region}++;
        }
        else {
            $uniques{$region} = 1;
        }
    }
    close DATA;
    return %uniques;
}


# setJointRegion($region_1_starting_position, $region_2_starting_position,
#   $window_size)
sub setJointRegion {
    my $self = shift;
    my $region_1_starting_position = shift;
    my $region_2_starting_position = shift;
    my $window_size = shift;
    my %unique_haplotypes;
    my $region;
    foreach my $line (keys %{$self->{"data"}}) {
        $region = substr(${$self->{'data'}}{$line}, $region_1_starting_position, $window_size) . 
            substr(${$self->{'data'}}{$line}, $region_2_starting_position, $window_size);
        if ($region !~ m/(0|1)+/) { # only allow 0's and 1's 
            next;
        }
        if (exists $unique_haplotypes{$region}) {
            $unique_haplotypes{$region}++;
        }
        else {
            $unique_haplotypes{$region} = 1;
        }
    }
    close DATA;
    return %unique_haplotypes;
}


1;
