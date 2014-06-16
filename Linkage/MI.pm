package Linkage::MI;
use Moose; # install the moose package

has 'locus_1' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'locus_2' => ( isa => 'Ref', is => 'rw',  required => 1);
has 'haplotype' => ( isa => 'Ref', is => 'rw',  required => 1);

sub log_2 {
	my $x = shift;
	return log($x)/log(2);
}

#currently returns output in nats (i.e. natural log rather than log base 2)
#accepts input 1 or 2
sub entropy {
    my $self = shift;
    my $locus = shift;
    my ($sum, $freq_1, $freq_2);
    if ($locus == 1) {
	$sum = $self->{locus_1}->{allele_a} + $self->{locus_1}->{allele_A};
	$freq_1 = $self->{locus_1}->{allele_a}/$sum;
	$freq_2 = $self->{locus_1}->{allele_A}/$sum;
    } 
    else {
	$sum = $self->{locus_2}->{allele_b} + $self->{locus_2}->{allele_B};
	$freq_1 = $self->{locus_2}->{allele_b}/$sum;
	$freq_2 = $self->{locus_2}->{allele_B}/$sum; 
    }
    return sprintf "%.4f", -1*($freq_1*log_2($freq_1) + $freq_2*log_2($freq_2));
}

#ab, aB, Ab, AB
sub joint_entropy {
    my $self = shift;
    my ($freq_1, $freq_2, $freq_3, $freq_4, $total_freq);
    $total_freq = $self->{haplotype}->{ab} + $self->{haplotype}->{aB} + $self->{haplotype}->{Ab} + $self->{haplotype}->{AB};
    $freq_1 = $self->{haplotype}->{ab}/$total_freq;
    $freq_2 = $self->{haplotype}->{aB}/$total_freq;
    $freq_3 = $self->{haplotype}->{Ab}/$total_freq;
    $freq_4 = $self->{haplotype}->{AB}/$total_freq;
    return sprintf "%.4f", -1 * ($freq_1*log_2($freq_1) + $freq_2*log_2($freq_2) + $freq_3*log_2($freq_3) + $freq_4*log_2($freq_4) );
}

sub mutual_information {
    my $self = shift;
    return sprintf "%.4f", ($self->entropy(1) + $self->entropy(2) - $self->joint_entropy() ); 
}

#return the MIR, or the normalized mutual information
sub mir {
    my $self = shift;
    my @sorted = sort {$a <=> $b} ($self->entropy(1), $self->entropy(2)); # sort numerically
    my $div = shift(@sorted); 
    return sprintf "%.4f", $self->mutual_information()/$div;
}

#attempt to calculate the error for two subregions, 
#with two haplotypes in each subregion
sub systematic_error {
    my $self = shift;
    my $m1 = 2;
    my $m2 = 2;
    my $m12 = $m1*$m2;
    my $N = $self->{haplotype}->{ab} + $self->{haplotype}->{aB} + $self->{haplotype}->{Ab} + $self->{haplotype}->{AB};
    return sprintf "%.4f", ($m12 - $m1 - $m2 + 1)/(2*$N);
}

#T statistic
#Null hypo: two regions are in linkage equilibrium
sub T_statistic {
    my $self = shift;
    my $N = $self->{haplotype}->{ab} + $self->{haplotype}->{aB} + $self->{haplotype}->{Ab} + $self->{haplotype}->{AB};
    return sprintf "%.4f", 2*$N*$self->mutual_information();
}

1;

=pod

=head1 Linkage::MI

Linkage::MI - Calculate mutual information of two polymorphic loci

=head1 SYNOPSIS

    use Linkage::MI;
 
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

