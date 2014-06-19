use Test::More tests => 2;
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
	reference => $reference
	);
my $bam_pileup = $bamutils->get_base_features(
	chr => $chr,
	pos => $pos
	);
$bamutils->print_base_features(
	pileup => $bam_pileup,
	output => $tempfile
	);
my $expected_file = "$Bin/example/example.pileup.txt";
compare_ok($tempfile, $expected_file, "pileup matches expected");
