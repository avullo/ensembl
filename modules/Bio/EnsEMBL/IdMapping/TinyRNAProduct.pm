=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::TinyRNAProduct - lightweight rnaproduct object

=head1 SYNOPSIS

  foreach my $rp (@{ $tr->get_all_RNAProducts() }) {
    my $lightweight_rp =
      Bio::EnsEMBL::IdMapping::TinyRNAProduct->new_fast( [
        $rp->dbID(),          $rp->stable_id(),
        $rp->version(),       $rp->created_date(),
        $rp->modified_date(), $tr->dbID(),
        $rp->seq(),
      ] );
  }

=head1 DESCRIPTION

This is a lightweight rnaproduct object for the stable Id mapping. See
the documentation in TinyFeature for general considerations about its
design.

=head1 METHODS

  transcript_id
  seq

=cut

package Bio::EnsEMBL::IdMapping::TinyRNAProduct;

# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  transcript_id
#  6  seq


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 transcript_id

  Arg[1]      : (optional) Int - the transcript internal Id ("dbID")
  Description : Getter/setter for the transcript internal Id this rnaproduct is
                attached to.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub transcript_id {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


=head2 seq

  Arg[1]      : (optional) String - the rnaproduct's sequence
  Description : Getter/setter for the rnaproduct's sequence.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}



1;

