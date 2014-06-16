package Linkage::multiMI;
use Moose; # install the moose package

has 'region_1' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'region_2' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'joint_haplotypes' => ( isa => 'Ref', is => 'rw',  required => 1);

sub log_2 {
	my $x = shift;
	return log($x)/log(2);
}

#currently returns output in log base 2
#accepts the locus number as the argument
sub entropy {
    my $self = shift;
    my $region_number = shift;
    my $region_tot_freq = 0;
    my %haplotypes = %{$self->{$region_number}};
    foreach my $key (keys %haplotypes) {
        $region_tot_freq += $haplotypes{$key};
    }
    my @proportions;
    foreach my $key (keys %haplotypes) {
        push (@proportions, $haplotypes{$key}/$region_tot_freq);
    }
    my $entro = 0;
    foreach my $prop (@proportions) {
        $entro += -1*($prop)*log_2($prop);
    }
    return sprintf "%.4f", $entro; }

#calculate the mutual information
sub mutual_information { 
    my $self = shift;
    my $m1 = scalar(keys %{$self->{"region_1"}}); #number of distinct haplotypes in region 1
    my $m2 = scalar(keys %{$self->{"region_2"}}); #number of distinct haplotypes in region 2
    my $m12 = $m1*$m2; #num distinct haplotypes in joint region
    my $region_1_tot_freq = 0;
    foreach my $key (keys %{$self->{"region_1"}}) {
        $region_1_tot_freq += ${$self->{"region_1"}}{$key};
    }
    my $region_2_tot_freq = 0;
    foreach my $key (keys %{$self->{"region_2"}}) {
        $region_2_tot_freq += ${$self->{"region_2"}}{$key};
    }
    my $total_joint_freq = 0;
    foreach my $key (keys %{$self->{"joint_haplotypes"}}) {
        $total_joint_freq += ${$self->{"joint_haplotypes"}}{$key};
    }
     
    my $sum_mi = 0;
    my $curr_mi = 0;
    my $joint_hap_prop;
    foreach my $r1hap (keys %{$self->{"region_1"}}) {
        foreach my $r2hap (keys %{$self->{"region_2"}}) {
            # print $r1hap, $r2hap, "\t", $self->{'joint_haplotypes'}{$r1hap . $r2hap}, "\n";
            if (exists $self->{'joint_haplotypes'}{$r1hap . $r2hap}) {
                $joint_hap_prop = $self->{'joint_haplotypes'}{$r1hap . $r2hap}/$total_joint_freq;
            }
            else {
                $joint_hap_prop = 0;
            }
            # loop through the joint haplotypes; 
            # each correct joint key is the concatenation of region1 and region2 keys
            if ($joint_hap_prop == 0) {
                $sum_mi += 0;
                # We define 0log0 = 0
            }
            else {
                $curr_mi = $joint_hap_prop * log_2( $joint_hap_prop / 
                    ( (${$self->{"region_1"}}{$r1hap}/$region_1_tot_freq)*(${$self->{"region_2"}}{$r2hap}/$region_2_tot_freq)  )
                    );
                $sum_mi += $curr_mi;
            }
        }
    }

    return "$sum_mi";
}

    

#return the MIR, or the normalized mutual information
sub mir {
    my $self = shift;
    my @entropies;
    push(@entropies, $self -> entropy("region_1"));
    push(@entropies, $self -> entropy("region_2"));
    my @sorted = sort {$a <=> $b} (@entropies); # sort numerically
    my $min = shift(@sorted); 
    return sprintf "%.4f", $self->mutual_information()/$min;
}

# Calculates the "systematic error" for two subregions 
# usage: systematic_error()
# note: not quite sure what N is supposed to be; guessed
sub systematic_error {
    my $self = shift;
    my $m1 = scalar(keys %{$self->{"region_1"}}); #number of distinct haplotypes in region 1
    my $m2 = scalar(keys %{$self->{"region_2"}}); #number of distinct haplotypes in region 2
    my $m12 = $m1*$m2; #num distinct haplotypes in joint region

    my $N = shift;
    return sprintf "%.4f", ($m12 - $m1 - $m2 + 1)/(2*$N);
}

# T statistic
# Null hypo: two regions are in linkage equilibrium
# note: not quite sure what N is
sub t_statistic {
    my $self = shift;
    my $N = shift;
    return sprintf "%.4f", 2*$N*$self->mutual_information();
}

1;

=pod

=head1 Linkage::multiMI

Linkage::MI - Calculate mutual information of two polymorphic loci

=head1 SYNOPSIS

    use Linkage::multiMI;
 
    my $data=Linkage::MI->new(locus_1=>{"allele_a" => 40, "allele_A" => 60}, locus_2=>{"allele_b" => 20, "allele_B" => 80}, haplotype=>{"ab"=>10, "aB"=>30, "Ab"=>10, "AB"=>50});
    my $h1=$data->entropy(1); # entropy at locus 1
    my $h2=$data->entropy(2); # entropy at locus 2
    my $mi=$data->entropy(); # mutual information


=head1 ABSTRACT

B<Linkage::MI> builds a object containing allele and haplotype frequencies and calculates entropy at each locus as well as mutual information of the two loci

=head1 DESCRIPTION

B<Linkage::MI> .....

=head1 METHODS

=head2 C<Constructor>

=over 4

=item
has 'locus_1' => ( isa => 'Ref', is => 'rw',  required => 1);

=item
has 'locus_2' => ( isa => 'Ref', is => 'rw',  required => 1);

=item
has 'haplotype' => ( isa => 'Ref', is => 'rw',  required => 1);

=back

=head2 C<entropy>

=over 4

=item Argument: 1 for locus 1; 2 for locus 2.

=item Return: Shannon Entropy

=back

=head1 DEPENDENCY

=over 4

=item Moose

=back

=head1 SEE ALSO

L<Linkage::xx>, L<Linkage::yy>

=head1 LICENSE

This software is released under the terms of the Academic Free License (AFL-3.0)
(http://www.opensource.org/licenses/AFL-3.0).

=head1 AUTHOR

Copyright (C) Dylan Sun
All rights Reserved

=cut

