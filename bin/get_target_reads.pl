#!/usr/bin/env perl

### get_lumpy_translocation_reads.pl ##############################################################
# Get the BAM entries containing potential translocation candidates from Lumpy output.

### HISTORY #######################################################################################
# Version       Date            Developer           Comments
# 0.01          2014-07-21      rdeborja            initial release

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Sam;
use File::Basename;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
	bam => undef,
	bedpe => undef
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
        "bedpe=s"
        ) or pod2usage(64);
    
    pod2usage(1) if $opts{'help'};
    pod2usage(-exitstatus => 0, -verbose => 2) if $opts{'man'};

    while(my ($arg, $value) = each(%opts)) {
        if (!defined($value)) {
            print "ERROR: Missing argument $arg\n";
            pod2usage(128);
            }
        }
    
    # read in all the targets from the BEDPE file (i.e. target file) into an array
    my @targets;
    open(my $target_fh, '<', $opts{'bedpe'});
    while(my $line = <$target_fh>) {
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;

        my @_targets = split(/\t/, $line);
        push(@targets, \@_targets);
    	}
    close($target_fh);

    # create a BAM file that will containg the filtered reads
    my $output_bam_filename = join('.',
        basename($opts{'bam'}, qw( .bam )),
        'targeted',
        'bam'
        );
    #my $output_bam = Bio::DB::Bam->open($output_bam_filename, "w");
    
    # create a file that will contain a list of read names
    # my $output_reads_list_filename = join('.',
    #     basename($opts{'bam'}, qw( .bam )),
    #     'readlist',
    #     'txt'
    #     );
    # open(my $ofh, '>', $output_reads_list_filename);

    # open the BAM file and loop through every entry and filtering for the targets from the
    # BEDPE data structure
    my $bam = Bio::DB::Bam->open($opts{'bam'});
    my $header = $bam->header();
    my $target_names = $header->target_name();

    # create a hash of read-ids, this will be used with Picard's FilterSamReads.jar to output
    # paired reads found in the 

    # write the header to the output file first
    # TO DO: add a PG line to include the current script used to filter/process the BAM file
    # my $header_write_status = $output_bam->header_write($header);

    # TODO: flip this around and make the main loop be the targets then search the BAM
    #       file for the current target
    # DONE!!!!! 
    # while(my $alignment = $bam->read1()) {
    #     my $chrA = $target_names->[$alignment->tid()];
    #     my $posA = $alignment->pos();
    #     my $chrB = $target_names->[$alignment->mtid()];
    #     my $posB = $alignment->mpos();
    #     my $read_name = $alignment->qname();

    #     # start looking for the positions in the translocation target array
    #     foreach my $trans_target (@targets) {
    #         if (($trans_target->[0] eq $chrA) & ($trans_target->[3] eq $chrB)) {
    #             if (($posA >= $trans_target->[1]) & ($posA <= $trans_target->[2])) {
    #                 if (($posB >= $trans_target->[4]) & ($posB <= $trans_target->[5])) {
    #                     #$output_bam->write1($alignment);
    #                     $_read_names{'$read_name'} = 1;
    #                     }
    #                 }
    #             }
    #         }
    #     }

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

    # foreach my $readname (keys(%_read_names)) {
    #     print {$ofh} "$readname\n";
    #     }
    # close($ofh);

    # currently manually running $PICARDROOT/FilterSamReads.jar to extract those reads that
    # contain the queryname in the read_list file

    return 0;
    }


__END__


=head1 NAME

get_lumpy_translocation_reads.pl

=head1 SYNOPSIS

B<get_lumpy_translocation_reads.pl> [options] [file ...]

    Options:
    --help          brief help message
    --man           full documentation
    --bam           BAM file to process
    --bedpe         regions containing potential translocations from lumpy

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--bam>

BAM file to process.

=item B<--bedpe>

BEDPE file containing regions for potential translocations.

=back

=head1 DESCRIPTION

B<get_lumpy_translocation_reads.pl> Get the BAM entries containing potential translocation candidates from Lumpy output.

=head1 EXAMPLE

get_lumpy_translocation_reads.pl

=head1 AUTHOR

Richard de Borja -- Molecular Genetics

The Hospital for Sick Children

=head1 SEE ALSO

=cut

