#!/usr/bin/perl

### piledriver.pl #################################################################################
# Print a table of features for a position for a BAM file.

### HISTORY #######################################################################################
# Version               Date                Developer       Comments
# 0.01                  2014-04-30          rdeborja        Initial development
# 0.02                  2014-06-16          rdeborja        updated for modified
#                                                           NGS::Tools::BamUtils library

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use NGS::Tools::BamUtils;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
    bam     => undef,
    chr     => undef,
    pos     => undef,
    ref     => '/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/hg19_random.fa',
    output  => ''
    );

### MAIN CALLER ###################################################################################
my $result = main();
exit($result);

### FUNCTIONS #####################################################################################

### main ##########################################################################################
# Description:
#   Main subroutine for program
# Input Variables:
#   %opts = command line arguments
# Output Variables:
#   N/A

sub main {
    # get the command line arguments
    GetOptions(
        \%opts,
        "help|?",
        'man',
        'bam|b=s',
        'chr|c=s',
        'pos|p=i',
        'ref|r=s',
        'output|o:s'
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
        bamfile     => $opts{'bam'},
        reference   => $opts{'ref'}
        );

    my $data = $bam->get_base_features(
        chr         => $opts{'chr'},
        pos         => $opts{'pos'}
        );
    $bam->print_base_features(pileup => $data, output => $args{'output'});

    return 0;
    }


__END__


=head1 NAME

piledriver.pl

=head1 SYNOPSIS

B<piledriver.pl> [options] [file ...]

    Options:
    --help          brief help message
    --man           full documentation
    --bam           BAM file to be processed
    --chr           chromosome/contig-name of interest
    --pos           position within reference genome of interest
    --ref           reference genome (default: hg19_random.fa)
    --output        filename to write data to, otherwise print to STDOUT (optional)

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--bam>

BAM file to be processed.

=item B<--chr>

Name of chromosome or contig as listed in the header of the BAM file.

=item B<--pos>

Position within the the reference genome to be processed.

=item B<--ref>

Reference genome (default: hg19_random.fa).

=item B<--output>

Filename to write data to.  Default is to write to STDOUT.

=back

=head1 DESCRIPTION

B<piledriver.pl> Print a table of features for a position for a BAM file.

=head1 EXAMPLE

piledriver.pl --bam file.bam --chr chr1 --pos 59570

=head1 AUTHOR

Richard de Borja

=head1 SEE ALSO

=cut

