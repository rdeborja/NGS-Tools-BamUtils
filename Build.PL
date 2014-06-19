use 5.006;
use strict;
use warnings FATAL => 'all';
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'NGS::Tools::BamUtils',
    license             => 'Artistic_2_0',
    dist_author         => q{Richard de Borja <richard.deborja@sickkids.ca>},
    dist_abstract       => 'Perl library for pileup using the Biosamtools package',
    dist_version_from   => 'lib/NGS/Tools/BamUtils.pm',
    release_status      => 'stable',
    configure_requires => {
        'Module::Build' => 0,
    },
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
	    'Bio::DB::Sam'   => 0,
			'Moose'          => 0
    },
    add_to_cleanup     => [ 'NGS-Tools-BamUtils-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
