package NGS::Tools::BamUtils::Roles::FilterReads;
use Moose::Role;
use MooseX::Params::Validate;

use strict;
use warnings FATAL => 'all';
use namespace::autoclean;
use autodie;
use Bio::DB::Sam;
use File::Temp qw(tempdir tempfile);

=head1 NAME

NGS::Tools::BamUtils::Roles::FilterReads

=head1 SYNOPSIS

A Perl Moose role for filtering reads based on a fixed set of criteria.

=head1 ATTRIBUTES AND DELEGATES

=head2 $obj->bam

BAM file to process

=cut

has 'bam' => (
    is          => 'rw',
    isa         => 'Str',
    reader		=> 'get_bam',
    writer		=> 'set_bam',
    required	=> 1
    );

=head1 SUBROUTINES/METHODS

=head2 $obj->get_discordant_read_pairs()

Retrieve the set of discordant pairs of reads from a BAM file.

=head3 Arguments:

=over 2

=item * bam: BAM file to process

=back

=cut

sub get_discordant_read_pairs {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		bam => {
			isa         => 'Str',
			required    => 0,
			default		=> $self->get_bam()
			},
		bedpe => {
			isa			=> 'Str',
			required	=> \1
			}
		);

	# get the targets
	my @targets;
	open(my $target_fh, '<', $args{'bedpe'});
	while(my $line = <$target_fh>) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;

		my @_targets = split(/\t/, $line);
		push(@targets, \@_targets);
		}
	close($target_fh);

	# open the BAM file for reading
	my $bam = Bio::DB::Bam->open($args{'bam'});

	# start looping through the targets and create a single BAM file for the predicted translocation of
	# interest
	foreach my $trans_target (@targets) {
		my %_readnames;
		while(my $alignment = $bam->read1()) {
	        my $chrA = $target_names->[$alignment->tid()];
	        my $posA = $alignment->pos();
	        my $chrB = $target_names->[$alignment->mtid()];
	        my $posB = $alignment->mpos();
	        my $readname = $alignment->qname();

            if (($trans_target->[0] eq $chrA) & ($trans_target->[3] eq $chrB)) {
                if (($posA >= $trans_target->[1]) & ($posA <= $trans_target->[2])) {
                    if (($posB >= $trans_target->[4]) & ($posB <= $trans_target->[5])) {
                        $_readnames{'$readname'} = 1;
                        }
                    }
                }
			}
		my $output_read_filename = join('.',
			basename($args{'bam'}, qw( .bam )),
			join('_',
				$trans_target->[0],
				$trans_target->[1],
				$trans_target->[2],
				$trans_target->[3],
				$trans_target->[4],
				$trans_target->[5]				
				),
			'readlist',
			'txt'
			);
		open(my $ofh, '>', $output_read_filename);
		foreach my $readname (keys(%_readnames)) {
			print {$ofh} "$readname\n";
			}
		close($ofh);
		}

	# my %return_values = (
	# 	output => $
	# 	);

	# return(\%return_values);
	return 0;
	}

=head2 $obj->get_translocation_reads()

Get discordant read pairs that align to different chromosomes.

=head3 Arguments:

=over 2

=item * bam: BAM file to process

=back

=cut

sub get_translocation_reads {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		bam => {
			isa         => 'Str',
			required    => 1
			},
		target => {
			isa			=> 'Str',
			required	=> 1
			}
		);

	my $bam = Bio::DB::Bam->open($args{'bam'});
	my $header = $bam->header();
	my $target_names = $header->target_name();
	while(my $alignment = $bam->read1()) {
		my $chr = $target_names->[$alignment->tid];
		my $chr_mate = $target_names->[$alignment->mtid];

		}


	my %return_values = (

		);

	return(\%return_values);
	}

=head1 AUTHOR

Richard de Borja, C<< <richard.deborja at sickkids.ca> >>

=head1 ACKNOWLEDGEMENT

Dr. Adam Shlien, PI -- The Hospital for Sick Children

Dr. Roland Arnold -- The Hospital for Sick Children

=head1 BUGS

Please report any bugs or feature requests to C<bug-test-test at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=test-test>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc NGS::Tools::BamUtils::Roles::FilterReads

You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=test-test>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/test-test>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/test-test>

=item * Search CPAN

L<http://search.cpan.org/dist/test-test/>

=back

=head1 ACKNOWLEDGEMENTS

=head1 LICENSE AND COPYRIGHT

Copyright 2013 Richard de Borja.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut

no Moose::Role;

1; # End of NGS::Tools::BamUtils::Roles::FilterReads
