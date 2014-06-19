package NGS::Tools::BamUtils::Roles::PileupPostProcess;
use Moose::Role;
use MooseX::Params::Validate;

use strict;
use warnings FATAL => 'all';
use namespace::autoclean;
use autodie;
use File::Slurp;
use Data::Dumper;

=head1 NAME

NGS::Tools::BamUtils::Roles::PileupPostProcess

=head1 SYNOPSIS

Post-process a pileup file.

=head1 ATTRIBUTES AND DELEGATES

=cut

my %_pileup_fields = (
	'Chr'				=> 0,
	'Pos'				=> 1,
	'Ref'				=> 2,
	'Call'				=> 3,
	'PhredScore'		=> 4,
	'ReadPos'			=> 5,
	'SAMFlag'			=> 6,
	'MapQual'			=> 7,
	'Strand'			=> 8,
	'NumMismatches'		=> 9,
	'MD'				=> 10,
	'ProperPair'		=> 11,
	'ISize'				=> 12,
	'ReadId'			=> 13,
	'overlap_reads'		=> 14,
	'overlap_consensus'	=> 15
	);

=head1 SUBROUTINES/METHODS

=head2 $obj->import_pileup_file()

Import the feature based pileup file.

=head3 Arguments:

=over 2

=item * pileup: argument

=back

=cut

sub import_pileup_file {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		pileup => {
			isa         => 'Str',
			required    => 1
			}
		);

	my $data = File::Slurp::read_file($args{'pileup'}, array_ref => 1);

	return($data)
	}

=head2 $obj->process_overlapping_reads()

For a given pileup file, process the overlapping reads.

=head3 Arguments:

=over 2

=item * data: pileup array reference generated from File::Slurp in import_pileup_file() method

=back

=cut

sub process_overlapping_reads {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		data => {
			isa         => 'ArrayRef',
			required    => 1
			}
		);

	my $data = $args{'data'};
	my @read_ids;

	# we will be recording duplicate read ids which will tell us if
	# read pairs overlap at the given position
	my %duplicate_read_ids;

	# need to create an array of read feature arrays to simplify the data structure
	my @_data;
	foreach my $read (@$data) {
		$read =~ s/^\s+//;
		$read =~ s/\s+$//;
		my @_tmp_data = split(/\t/, $read);
		push(@_tmp_data, qw(FALSE NA));
		push(@_data, [@_tmp_data]);
		# $duplicate_read_ids{$_tmp_data[$_pileup_fields{'ReadId'}]}++;
		# $duplicate_read_ids{$_tmp_data[$_pileup_fields{'ReadId'}]} = {'base' => []};
		push(@{ $duplicate_read_ids{$_tmp_data[$_pileup_fields{'ReadId'}]}->{'base'} }, $_tmp_data[$_pileup_fields{'Call'}]);
		$duplicate_read_ids{$_tmp_data[$_pileup_fields{'ReadId'}]}->{'count'}++;
		}

	# remove any key/value pairs in duplicate_read_ids if the value is not 2
	while(my ($key, $value) = each(%duplicate_read_ids)) {
		# if ($value < 2) {
		if ($value->{'count'} < 2) {
			delete($duplicate_read_ids{$key});
			}
		}
	# determine whether Read Ids are in duplicate meaning pairs of reads that overlap
	for(my $i = 0; $i < scalar(@_data); $i++) {
		if (exists($duplicate_read_ids{$_data[$i][$_pileup_fields{'ReadId'}]})) {
			$_data[$i][$_pileup_fields{'overlap_reads'}] = 'TRUE';
			my %_base_hash;
			# for those reads that are in duplicate, check that the "Call" bases are the same, if they are then
			# let overlap_consensus be TRUE otherwise FALSE, all non-overlapping reads will remain NA
			foreach my $base (@{ $duplicate_read_ids{$_data[$i][$_pileup_fields{'ReadId'}]}->{'base'} }) {
				unless ($_base_hash{$base}++) {
					$_data[$i][$_pileup_fields{'overlap_consensus'}] = 'FALSE'
					}
				else {
					$_data[$i][$_pileup_fields{'overlap_consensus'}] = 'TRUE'
					}
				}
			}
		}

	return(\@_data);
	}

=head2 $obj->sort_array_by_readid()

A method to sort an array of arrays by the readid column.  Reference: 
http://www.sitepoint.com/forums/showthread.php?659680-Perl-Sort-Array-of-Arrays

=head3 Arguments:

=over 2

=item * data: array reference to sort by ReadId

=back

=head3 Return Values:

=over 2

=item * \@_sorted_data: reference to the sorted array

=back

=cut

sub sort_array_by_readid {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		data => {
			isa         => 'ArrayRef',
			required    => 1
			}
		);

	my $data = $args{'data'};
	my @_sorted_data = sort {$a->[$_pileup_fields{'ReadId'}] cmp $b->[$_pileup_fields{'ReadId'}]} @{$data};

	return(\@_sorted_data);
	}

=head2 $obj->header()

Get the pileup header from the pileup file and add overlap_reads.

=head3 Arguments:

=over 2

=item * data: pileup data structure

=back

=head3 Return Values:

=over 2

=item * \@header: array reference containing the pileup file header

=back

=cut

sub header {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		data => {
			isa         => 'ArrayRef',
			required    => 1
			}
		);
	my $data = $args{'data'};
	my $header_line = shift(@$data);
	$header_line =~ s/^\s+//;
	$header_line =~ s/\s+$//;
	my @header = split(/\t/, $header_line);
	push(@header, qw(overlap_reads overlap_consensus));

	return(\@header);
	}

=head2 $obj->print_pileup_data()

Write the pileup data to a file.

=head3 Arguments:

=over 2

=item * data: data structure containing pileup information.

=back

=cut

sub print_pileup_data {
	my $self = shift;
	my %args = validated_hash(
		\@_,
		data => {
			isa         => 'ArrayRef',
			required    => 1
			},
		output => {
			isa			=> 'Str',
			required	=> 0,
			default		=> ''
			},
		header => {
			isa			=> 'ArrayRef',
			required	=> 1
			}
		);

	my $header = $args{'header'};

	# print to the file only if an output file was provided, otherwise print to STDOUT.
	my $ofh;
	if ($args{'output'} ne '') {
		open($ofh, '>', $args{'output'});
		select($ofh);
		}
	print join("\t", @{$header}), "\n";
	foreach my $read (@{$args{'data'}}) {
		print join("\t",
			$read->[$_pileup_fields{'Chr'}],
			$read->[$_pileup_fields{'Pos'}],
			$read->[$_pileup_fields{'Ref'}],
			$read->[$_pileup_fields{'Call'}],
			$read->[$_pileup_fields{'PhredScore'}],
			$read->[$_pileup_fields{'ReadPos'}],
			$read->[$_pileup_fields{'SAMFlag'}],
			$read->[$_pileup_fields{'MapQual'}],
			$read->[$_pileup_fields{'Strand'}],
			$read->[$_pileup_fields{'NumMismatches'}],
			$read->[$_pileup_fields{'MD'}],
			$read->[$_pileup_fields{'ProperPair'}],
			$read->[$_pileup_fields{'ISize'}],
			$read->[$_pileup_fields{'ReadId'}],
			$read->[$_pileup_fields{'overlap_reads'}],
			$read->[$_pileup_fields{'overlap_consensus'}],
			), "\n";
		}
	if (defined($ofh)) {
		select(STDOUT);
		close($ofh);
		}

	return 0;
	}

=head1 AUTHOR

Richard de Borja, C<< <richard.deborja at sickkids.ca> >>

=head1 ACKNOWLEDGEMENT

Dr. Adam Shlien, PI -- The Hospital for Sick Children

Dr. Roland Arnold -- The Hospital for Sick Children

Dr. M. Anaka -- The Hospital for Sick Children

Andrej Rosic -- The Hospital for Sick Children

=head1 BUGS

Please report any bugs or feature requests to C<bug-test-test at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=test-test>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc NGS::Tools::BamUtils::Roles::PileupPostProcess

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

1; # End of NGS::Tools::BamUtils::Roles::PileupPostProcess
