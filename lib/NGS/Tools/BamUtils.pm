package NGS::Tools::BamUtils;
use Moose;
use MooseX::Params::Validate;

with 'NGS::Tools::BamUtils::Roles::PileupUtilities';
with 'NGS::Tools::BamUtils::Roles::PileupPostProcess';

use 5.006;
use strict;
use warnings;
use autodie;
use Bio::DB::Sam;
use File::Temp qw(tempdir tempfile);

=head1 NAME

NGS::Tools::BamUtils

A suite of utilities to perform various tasks on an indexed BAM file.

=head1 VERSION

Version 0.06

=cut

our $VERSION = '0.06';


=head1 SYNOPSIS

The BamUtils module extends the BioSamtools framework to perform common
and simple tasks.

Usage:

    use NGS::Tools::BamUtils;

    my $obj = NGS::Tools::BamUtils->new();
    $obj->get_pileup(chr => $chr, pos => $pos);
    $obj->get_nucleotide_counts(chr => $chr, pos => $pos);
    ...

=head1 INSTANCE VARIABLES

=head1 SUBROUTINES/METHODS

=item * $args: reference to hash of arguments

=back

=cut

sub BUILD {
    my $self = shift;
    my $args = shift;
    }

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc NGS::Tools::BamUtils


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=NGS-Tools-BamUtils>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/NGS-Tools-BamUtils>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/NGS-Tools-BamUtils>

=item * Search CPAN

L<http://search.cpan.org/dist/NGS-Tools-BamUtils/>

=back


=head1 ACKNOWLEDGEMENTS

Michelle Chan-Seng-Yue C<< <michelle.chansengyue at oicr.on.ca > >>

=head1 LICENSE AND COPYRIGHT

Copyright 2014 Richard de Borja.

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

no Moose;

__PACKAGE__->meta->make_immutable;

1; # End of NGS::Tools::BamUtils
