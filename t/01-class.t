use Test::More tests => 1;
use Test::Files;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Temp qw(tempfile tempdir);
use Data::Dumper;

BEGIN: {
    use_ok('NGS::Tools::BamUtils');
	}

my $bamfile = "$Bin/example/exampleBAM.bam";
my $reference = "$Bin/example/exampleFASTA.fasta";
my $chr = 'chr1';
my $pos = 255;
my (undef, $tempfile) = tempfile(DIR => '.', UNLINK => 1);
my $bamutils = NGS::Tools::BamUtils->new(
	bamfile => $bamfile,
	reference => $reference,
	output => $tempfile
	);
my @bam_pileup = $bamutils->get_pileup(
	chr => $chr,
	pos => $pos
	);
$bamutils->print_base_frequency_table(
	pileup => \@bam_pileup
	);