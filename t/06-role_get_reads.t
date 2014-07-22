use Test::More tests => 1;
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
my $test_class = $test_class_factory->class_for('NGS::Tools::BamUtils::Roles::FilterReads');

# instantiate the test class based on the given role
my $bamutils;
my $bam = "$Bin/example/exampleBam.bam";

lives_ok
	{
		$bamutils = $test_class->new(
			bam => $bam
			);
		}
	'Class instantiated';

$reads = $bamutils->get_translocation_reads();
