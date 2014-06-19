package NGS::Tools::BamUtils::Roles::PileupUtilities;
use Moose::Role;
use MooseX::Params::Validate;

use strict;
use warnings FATAL => 'all';
use namespace::autoclean;
use autodie;
use Bio::DB::Sam;

=head1 NAME

NGS::Tools::BamUtils::Roles::PileupUtilities

=head1 SYNOPSIS

A Perl Moose Role for handling the creation of a modified pileup file.

=head1 ATTRIBUTES AND DELEGATES

=head2 $obj->bamfile

BAM file to be processed, defined in object constructor.

=cut

has 'bamfile' => (
    is          => 'rw',
    isa         => 'Str',
    reader      => 'get_bamfile',
    writer      => 'set_bamfile'
    );

=head2 $obj->reference

Reference genome (mandatory field)

=cut

has 'reference' => (
    is          => 'rw',
    isa         => 'Str',
    reader      => 'get_reference',
    writer      => 'set_reference'
    );

# setup a base frequency table
our %_base_frequency = (
    'A'         => 0,
    'T'         => 0,
    'C'         => 0,
    'G'         => 0,
    'N'         => 0
    );

our %_truth_table = (
    0           => 'FALSE',
    1           => 'TRUE'
    );

=head1 SUBROUTINES/METHODS

=head2 get_pileup(chr => $chr, pos => $pos, qual => $base_qual)

=head3 Description

Get the pileup of bases at a single position in a BAM file. Default qual = 0

=head3 Arguments:

=over 2

=item * chr: chromosome/contig for position of interest

=item * pos: position within chromosome/contig of interest

=item * qual: base quality (Phred scale) for read filtering

=item * max_reads: maximum number of reads to collect for pileup

=back

=head3 Return Values

=over 2

=item * base_pileup: returns an array containing the base pileup

=back

=cut

sub get_pileup {
    my $self = shift;
    my %args = validated_hash(
        \@_,
        bam => {
            isa         => 'Str',
            required    => 0,
            default     => $self->get_bamfile()
            },
        reference => {
            isa         => 'Str',
            required    => 0,
            default     => $self->get_reference()
            },
        chr => {
            isa         => 'Str',
            required    => 1
            },
        pos => {
            isa         => 'Int',
            required    => 1
            },
        qual => {
            isa         => 'Int',
            required    => 0,
            default     => 0
            },
        max_reads => {
            isa         => 'Int',
            required    => 0,
            default     => 100000000
            }
        );

    # create a new BioSamtools object
    my $bam = Bio::DB::Sam->new(
        -bam        => $args{'bam'},
        -fasta      => $args{'reference'},
        );
    $bam->max_pileup_cnt($args{'max_reads'});
    my %base_frequency = %_base_frequency;
    my $region = "$args{'chr'}\:$args{'pos'}-$args{'pos'}";
    my @base_pileup;
    my $callback = sub {
        my ($seqid, $pos, $pileup_array) = @_;
        for my $pileup (@$pileup_array) {
            my $b = $pileup->alignment;
            next if $pileup->indel or $pileup->is_refskip;
            my $qbase = substr($b->qseq, $pileup->qpos, 1);
            next if $qbase =~ m/[nN]/;
            my $qscore = $b->qscore->[$pileup->qpos];
            next unless $qscore > $args{'qual'};
            if ($args{'pos'} eq $pos) {
                push(@base_pileup, $qbase);
                }
            }

        };
    $bam->pileup($region, $callback);
    return(@base_pileup);
    }

=head2 $obj->get_base_features(chr => $chr, pos => $pos)

Get a set of features for a given chromosome and position.

=head3 Arguments:

=over 2

=item * chr: chromosome/contig for position of interest

=item * pos: position within chromosome/contig of interest

=item * reference: reference genome in FASTA format

=item * notes: 

=back

=cut

sub get_base_features {
    my $self = shift;
    my %args = validated_hash(
        \@_,
        bam => {
            isa         => 'Str',
            required    => 0,
            default     => $self->get_bamfile
            },
        chr => {
            isa         => 'Str',
            required    => 1
            },
        pos => {
            isa         => 'Int',
            required    => 1
            },
        reference => {
            isa         => 'Str',
            required    => 0,
            default     => $self->get_reference()
            },
        qual => {
            isa         => 'Int',
            required    => 0,
            default     => 0
            },
        append => {
            isa         => 'Int',
            required    => 0,
            default     => 0
            },
        max_reads => {
            isa         => 'Int',
            required    => 0,
            default     => 100000000
            }
        );

    # store each pileup line as an entry in the array, we will return this
    # array as a reference
    my @pileup_array;
    my $bam;
    if (defined($args{'reference'})) {
        $bam = Bio::DB::Sam->new(
            -bam            => $args{'bam'},
            -fasta          => $args{'reference'},
            -expand_flags   => 1
            );
    }
    else {
        $bam = Bio::DB::Sam->new(
            -bam            => $args{'bam'},
            -expand_flags   => 1
            );
        }
    $bam->max_pileup_cnt($args{'max_reads'});
    my $region = "$args{'chr'}\:$args{'pos'}-$args{'pos'}";
    my @base_features;
    my @base_pileup;
    my @quality_pileup;
    my $callback = sub {
        my ($seqid, $pos, $pileup_array) = @_;
        my $refbase = $bam->segment($seqid,$pos,$pos)->dna;
        for my $pileup (@$pileup_array) {
            my $b = $pileup->alignment;
            my $strand;
            if ($b->strand == 1) {
                $strand = '+';
                }
            elsif ($b->strand == -1) {
                $strand = '-';
                }
            # the tags in a SAM file are obtained as a tab separated string
            my @tags = split(/\t/, $b->aux);
            my $mismatches = 'NA';
            my $mdtag = 'NA';
            foreach my $tag (@tags) {
                if ($tag =~ m/^NM/) {
                    # an NM tag looks like the following: NM:i:0
                    (undef, undef, $mismatches) = split(/\:/, $tag);
                    }
                elsif ($tag =~ m/^MD/) {
                    # a MD tag looks like the following: MD:Z:151M
                    (undef, undef, $mdtag) = split(/\:/, $tag);
                    }
                }
            my $readid = $b->qname;
            next if $pileup->indel or $pileup->is_refskip;
            my $qbase = substr($b->qseq, $pileup->qpos, 1);
            next if $qbase =~ m/[nN]/;
            my $qscore = $b->qscore->[$pileup->qpos];
            next unless $qscore >= $args{'qual'};
            if ($args{'pos'} eq $pos) {
                my %_pileup_features = (
                    'base'          => $qbase,
                    'qual'          => $qscore,
                    'readpos'       => $pileup->pos,
                    'alignpos'      => $pos
                    );                
                push(@base_pileup, $qbase);
                push(@quality_pileup, $qscore);
                $_base_frequency{$qbase}++;

                # print join("\t",
                my $pileup_entry = join("\t",
                    $seqid,
                    $pos,
                    $refbase,
                    $qbase,
                    $qscore,
                    $pileup->pos,
                    $b->flag(),
                    $b->qual(),
                    $strand,
                    $mismatches,
                    $mdtag,
                    $_truth_table{$b->proper_pair()},
                    $b->isize(),
                    $readid
                    );
#                    ) . "\n";
                push (@pileup_array, $pileup_entry);
                }
            }

        };

    # add a file header
    if ($args{'append'} == 0) {
        # $self->_print_base_features_header();
        push(@pileup_array, $self->_get_base_features_header());
        }
    else {
        # don't add a header
        }
    $bam->pileup($region, $callback);

    return(\@pileup_array);
    # return(0);
    }

=head2 $obj->_print_base_features_header()

Print the column header for the get_base_features() output.

=head3 Arguments:

=over 2

=item * notes: print the note column header if requested

=back

=cut

sub _print_base_features_header {
    my $self = shift;
    print join("\t",
        'Chr',
        'Pos',
        'Ref',
        'Call',
        'PhredScore',
        'ReadPos',
        'SAMFlag',
        'MapQual',
        'Strand',
        'NumMismatches',
        'MD',
        'ProperPair',
        'ISize',
        'ReadId'
        ) . "\n";
    }

=head2 $obj->print_base_features()

Print the base features array reference generated from get_base_features().

=head3 Arguments:

=over 2

=item * pileup: an array reference containing the output from get_base_features().

=back

=cut

sub print_base_features {
    my $self = shift;
    my %args = validated_hash(
        \@_,
        pileup => {
            isa         => 'ArrayRef',
            required    => 1
            },
        output => {
            isa         => 'Str',
            required    => 0,
            default     => ''
            }
        );

    # open a file for outputting only if $args{'output'} is passed to the method,
    # otherwise print to STDOUT
    my $ofh;
    if ($args{'output'} ne '') {
        open($ofh, '>', $args{'output'});
        select($ofh);
        }

    foreach my $pileup_line (@{$args{'pileup'}}) {
        print $pileup_line, "\n";
        }

    # if output filename provided, close the file and return select to STDOUT
    if ($args{'output'} ne '') {
        close($ofh);
        select(STDOUT);
        }

    return(0);
    }


=head2 $obj->_get_base_features_header()

Return the header line for the base features.

=head3 Return Values:

=over 2

=item * header: a tab separated string containing the base header line

=back

=cut

sub _get_base_features_header {
    my $self = shift;

    return(
        join("\t",
            'Chr',
            'Pos',
            'Ref',
            'Call',
            'PhredScore',
            'ReadPos',
            'SAMFlag',
            'MapQual',
            'Strand',
            'NumMismatches',
            'MD',
            'ProperPair',
            'ISize',
            'ReadId'
            )
        );
    }

=head2 $obj->_print_nucleotide_counts_header()

Print the header line for the nucleotide counts table.

=cut

sub _print_nucleotide_counts_header {
    my $self = shift;
    print join("\t",
        'A',
        'T',
        'C',
        'G',
        'N'
        ) . "\n";
    return 0;
    }

=head2 get_nucleotide_counts(pileup => $pileup_array_ref)

Get the nucleotide counts for the given position.  Returns a hash ref
with the summary of base counts.

=head3 Arguments

=over 2

=item * pileup: Array containing nucleotide bases generated from get_pileup().

=back

=cut

sub get_nucleotide_counts {
    my $self = shift;
    my %args = validated_hash(
    	\@_,
        pileup => {
        	isa			=> 'ArrayRef',
        	required	=> 1
        	}
        );
    my %base_frequency = %_base_frequency;
    foreach my $query_base (@{$args{'pileup'}}) {
        $base_frequency{$query_base}++;
        }

    return(\%base_frequency);
    }

=head2 print_nucleotide_counts(base_frequency => $hash_ref)

Print the nucleotide counts generated by get_nucleotide_counts method

=head3 Arguments:

=over 2

=item * base_frequency: a hash reference containing the base frequency data from the get_nucleotide_counts method

=back

=cut

sub print_nucleotide_counts {
    my $self = shift;
    my %args = validated_hash(
    	\@_,
        base_frequency => {
        	isa			=> 'HashRef',
        	required	=> 1
        	}
        );

    # use the class variable declared above
    my $base_frequency = $args{'base_frequency'};

    print join("\t",
        $base_frequency->{'A'},
        $base_frequency->{'T'},
        $base_frequency->{'C'},
        $base_frequency->{'G'},
        $base_frequency->{'N'}
        ) . "\n";

    return 0;
    }

=head2 $obj->print_base_frequency_table(pileup => $pileup_array_ref)

A macro to print the base frequency table based on the pileup for a given position.

=head3 Arguments

=over 2

=item * pileup: Array containing nucleotide bases generated from get_pileup().

=back

=cut

sub print_base_frequency_table {
    my $self = shift;
    my %args = validated_hash(
        \@_,
        pileup => {
        	isa			=> 'ArrayRef',
        	required	=> 1
        	}
        );
    my $pileup = $args{'pileup'};
    $self->_print_nucleotide_counts_header();
    my $base_frequency = $self->get_nucleotide_counts(pileup => $pileup);
    $self->print_nucleotide_counts(base_frequency  => $base_frequency);

    return 0;
    }


=head1 AUTHOR

Richard de Borja, C<< <richard.deborja at sickkids.ca> >>

=head1 ACKNOWLEDGEMENT

Dr. Adam Shlien, PI -- The Hospital for Sick Children

Dr. Roland Arnold -- The Hospital for Sick Children

Andrej Rosic -- The Hospital for Sick Children / Waterloo University

=head1 BUGS

Please report any bugs or feature requests to C<bug-test-test at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=test-test>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc NGS::Tools::BamUtils::Roles::PileupUtilities

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

1; # End of NGS::Tools::BamUtils::Roles::PileupUtilities
