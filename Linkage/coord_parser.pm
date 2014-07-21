package Linkage::coord_parser;
use Moose; 

has 'location_file' => ( isa => 'Any', is => 'rw',  required => 1);
has 'sites_file' => ( isa => 'Any', is => 'rw',  required => 1);
has 'locations' => ( isa => 'ArrayRef', is => 'rw',  required => 0);
has 'site_info' => ( isa => 'ArrayRef', is => 'rw',  required => 0);
has 'location_info' => ( isa => 'ArrayRef', is => 'rw',  required => 0);
has 'data' => ( isa => 'Any', is => 'rw',  required => 0);

sub parseLocations {
    my $self = shift;
    open(DATA, $self->{"location_file"}) || die "Can't open snp file: $!";
    my $first_line = <DATA>;
    @{$self->{"location_info"}} = split (/\s+/, $first_line);
    while (my $line = <DATA>) {
        chomp $line;
        push(@{$self->{"locations"}}, $line);
    }
    close DATA;
}

sub parseSites {
    my $self = shift;
    my $counter = 0;
    my $string = "";
    my $is_info = 0;
    open(DATA, $self->{"sites_file"}) || die "Can't open snp file: $!";
    my $first_line = <DATA>;
    @{$self->{"site_info"}} = split (/\s+/, $first_line);
    while (my $line = <DATA>) {
        chomp $line;
        if ($line =~ m/^S/) {
            ${$self->{"data"}}{$counter} = $string;
            $counter++;
            $string = "";
        }
        else {
            $string = $string . $line;
        }
    }
    close DATA;
    return $counter; 

}

# getRegion($starting_position, $window_size)
sub getRegion {
    my $self = shift;
    my $starting_position = shift;
    my $window_size = shift;
    my %uniques;
    my $region;
    my $first_snp = -1;
    my $num_sites = 0;
    my $counter = 0;
    # for (my $i; $i < ${$self -> {"location_info"}}[1]; $i++)
    foreach my $snp_loc (@{$self->{"locations"}}) {
        if ($starting_position <= $snp_loc && $snp_loc <= $starting_position + $window_size - 1) {
            if ($first_snp < 0 ) {
                $first_snp = $counter;
            }
            $num_sites++;
        }
        elsif ($snp_loc > $starting_position + $window_size - 1) {
           last;
       }
       $counter++;
    }
    foreach my $line (keys %{$self->{"data"}}) {
        $region = substr(${$self->{'data'}}{$line}, $first_snp, $num_sites);
        # so at every window with only 1 snp site, there are only 2 haplotypes?
        if (exists $uniques{$region}) {
            $uniques{$region}++;
        }
        else {
            $uniques{$region} = 1;
        }
    }
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
