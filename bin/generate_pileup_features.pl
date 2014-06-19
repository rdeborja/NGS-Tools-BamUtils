#!/usr/bin/perl

### generate_pileup_features.pl ####################################################################
# Generate a file containing pileup data and BAM file features.

### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2014-06-14      rdeborja            initial development

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use NGS::Tools::BamUtils;
use File::Temp qw(tempdir tempfile);

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
	bam => undef,
	reference => undef,
	output => '',
    chr => undef,
    pos => undef
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
        "man",
        "bam|b=s",
        "reference|r=s",
        "output|o:s",
        "chr|c=s",
        "pos|p=s"
        ) or pod2usage(64);
    
    pod2usage(1) if $opts{'help'};
    pod2usage(-exitstatus => 0, -verbose => 2) if $opts{'man'};

    while(my ($arg, $value) = each(%opts)) {
        if (!defined($value)) {
            print "ERROR: Missing argument \n";
            pod2usage(128);
            }
        }

    my $bam = NGS::Tools::BamUtils->new(
        bamfile     => $opts{'bam'},
        reference   => $opts{'reference'}
        );

    # create the pileup array
    my $data = $bam->get_base_features(
        chr         => $opts{'chr'},
        pos         => $opts{'pos'}
        );

    # perform the overlap analysis to include overlap_reads and overlap_consensus
    # covariates to the pileup file
    my $header = $bam->header(
    	data => $data
    	);
    my $pileup_data = $bam->process_overlapping_reads(
    	data => $data
    	);
    my $sorted_pileup_data = $bam->sort_array_by_readid(
    	data => $pileup_data
    	);
    $bam->print_pileup_data(
    	data => $sorted_pileup_data,
    	header => $header,
    	output => $opts{'output'}
    	);

    return 0;
    }


__END__


=head1 NAME

generate_pileup_features.pl

=head1 SYNOPSIS

B<generate_pileup_features.pl> [options] [file ...]

    Options:
    --help          brief help message
    --man           full documentation
    --bam           BAM file to process (required)
    --reference     FASTA reference genome (required)
    --output        filename to write data to (optional)
    --chr           chromsome/contig to process
    --pos           position within chromosome/contig to process

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--bam>

BAM file to process.

=item B<--reference>

FASTA reference genome used with BAM file.

=item B<--output>

Filename to write data to.  If left blank, output will default to STDOUT.

=item B<--chr>

Name of chromosome/contig to query.

=item B<--pos>

Position within chromosome/contig defined in --chr to query.

=back

=head1 DESCRIPTION

B<generate_pileup_features.pl> Generate a file containing pileup data and BAM file features.

=head1 EXAMPLE

generate_pileup_features.pl --bam file.bam --reference ref.fasta --output output.txt --chr chr1 --pos 100

=head1 AUTHOR

Richard de Borja -- Molecular Genetics

The Hospital for Sick Children

=head1 SEE ALSO

=cut

