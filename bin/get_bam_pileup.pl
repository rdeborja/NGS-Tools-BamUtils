#!/usr/bin/perl

### get_bam_pileup.pl #############################################################################
# Get the base pileup for a single position in a BAM file.  The BAM file must be indexed.

### HISTORY #######################################################################################
# Version           Date            Developer           Comments
# 0.01              2014-04-30      rdeborja            Initial version.

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);
use Cwd;
use NGS::Tools::BamUtils;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
    bam         => undef,
    chr         => undef,
    pos         => undef,
    qual        => 0
    );

### MAIN CALLER ###################################################################################
my $result = main();
exit($result);

### FUNCTIONS #####################################################################################

### main ##########################################################################################
# Description:
#    Main subroutine for program
# Input Variables:
#    %opts = command line arguments
# Output Variables:
#    N/A

sub main {
    # get the command line arguments
    GetOptions(
        \%opts,
        "help|?",
        "man",
        "bam|b=s",
        "chr|c=s",
        "pos|p=s",
        "qual|q=s",
        ) or pod2usage(64);

    if ($opts{'help'}) { pod2usage(1) };
    if ($opts{'man'}) { pod2usage(-exitstatus => 0, -verbose => 2) };

    while(my ($arg, $value) = each(%opts)) {
        if (!defined $value) {
            print "ERROR: Missing argument $arg\n";
            pod2usage(128);
            }
        }

    my $bam = NGS::Tools::BamUtils->new(
        bamfile     => $opts{'bam'}
        );

    my @bam_pileup = $bam->get_pileup(
        chr => $opts{'chr'},
        pos => $opts{'pos'},
        qual => $opts{'qual'}
        );
    $bam->print_base_frequency_table(
        pileup => \@bam_pileup
        );
    return 0;

    }


__END__


=head1 NAME



=head1 SYNOPSIS

B<> [options] [file ...]

    Options:
    --help          brief help message
    --man           full documentation
    --bam           BAM file to process
    --chr           chromosome/contig for pileup
    --pos           position within chromosome/contig to generate pileup
    --qual          base quality minimum value for inclusion into counts (default: 0)

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--bam>

BAM file to be used for processing.

=item B<--chr>

Chromosome/contig to be used within BAM file.

=item B<--pos>

Position within in the chromosome/contig to be processed.

=item B<--qual>

Minimum base quality for base within read to be processed. Default 0

=back

=head1 DESCRIPTION

B<get_bam_pileup.pl> Get the BAM pileup for a specific genomic position
and identify very base for a read.

=head1 EXAMPLE

get_bam_pileup --bam example.bam --chr chr1 --pos 59570 --qual 20

=head1 AUTHOR

Richard de Borja

=head1 SEE ALSO

=cut


