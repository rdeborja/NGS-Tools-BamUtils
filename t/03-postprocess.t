use Test::More tests => 2;
use Test::Moose;
use Test::Exception;
use MooseX::ClassCompositor;
use Test::Files;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Temp qw(tempfile tempdir);
use Data::Dumper;

# setup the class creation process
my $test_class_factory = MooseX::ClassCompositor->new(
	{ class_basename => 'Test' }
	);

# create a temporary class based on the given Moose::Role package
my $test_class = $test_class_factory->class_for('NGS::Tools::BamUtils::Roles::PileupPostProcess');

# instantiate the test class based on the given role
my $pileup;
lives_ok
	{
		$pileup = $test_class->new();
		}
	'Class instantiated';


# my $file = "$Bin/example/pileup.txt";
my $file = "$Bin/example/example.pileup.txt";

my $data = $pileup->import_pileup_file(pileup => $file);
my $header = $pileup->header(data => $data);

my $pileup_data = $pileup->process_overlapping_reads(
	data => $data
	);
$sorted_pileup_data = $pileup->sort_array_by_readid(
	data => $pileup_data
	);
my $tempdir = File::Temp::tempdir(DIR => '.', CLEANUP => 1);
my (undef, $tempfile) = File::Temp::tempfile(DIR => $tempdir);
$pileup->print_pileup_data(
	data => $sorted_pileup_data,
	header => $header,
	output => $tempfile
	);
my $expected_file = "$Bin/example/03_expected_features.txt";
compare_ok($tempfile, $expected_file, "pileup features matches expected");