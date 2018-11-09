=head1 LICENSE

See the NOTICE file distributed with this work for additional
information regarding copyright ownership.

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



package XrefParser::UniProtParser::Transformer;

use strict;
use warnings;

use Carp;
use Readonly;

# FIXME: for testing
use Data::Dumper;


# FIXME: this should probably be combined with
# Extractor::%supported_taxon_database_qualifiers to make sure
# database qualifiers stay in sync
Readonly my %taxonomy_ids_from_taxdb_codes
  => {
      # NCBI taxon codes and Ensembl taxonomy IDs are identical
      'NCBI_TaxID' => sub { return $_[0]; },
    };

Readonly my %whitelisted_crossreference_sources
  => (
      'ChEMBL'  => 1,
      'EMBL'    => 1,
      'Ensembl' => 1,
      'MEROPS'  => 1,
      'PDB'     => 1,
    );

Readonly my $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD => 2;
Readonly my %source_selection_criteria_for_status
  => (
      'Reviewed'   => [ 'Uniprot/SWISSPROT',
                        sub {
                          return 'sequence_mapped';
                        }, ],
      'Unreviewed' => [ 'Uniprot/SPTREMBL', sub {
                          my ( $level ) = @_;
                          return ( $level <= $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD ) ?
                            'sequence_mapped' :
                            "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD";
                        }, ],
    );



sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
              'dbh'        => $arg_ref->{'dbh'},
              'species_id' => $arg_ref->{'species_id'},
              'maps'       => {},
            };
  my $class = ref $proto || $proto;
  bless $self, $class;

  $self->_load_maps();

  return $self;
}


# Transforms extracted record into form that can be consumed by
# BaseParser::upload_xref_object_graphs().
sub transform {
  my ( $self, $extracted_record ) = @_;

  $self->{'extracted_record'} = $extracted_record;

  # Only proceed if at least one taxon code in the entry maps
  # to a valid Ensembl species ID
  my $xref_multiplicity = $self->_recognised_taxon_ids();
  if ( ! $xref_multiplicity ) {
    return;
  }

  # Skip unreviewed entries
  if ( $self->_entry_is_unreviewed() ) {
    return;
  }

  my ( $accession, @synonyms )
    = @{ $extracted_record->{'accession_numbers'} };

  my $xref_graph_node
    = {
       'ACCESSION'     => $accession,
       'DESCRIPTION'   => $extracted_record->{'description'},
       'INFO_TYPE'     => 'SEQUENCE_MATCH',  # FIXME: is this always right?
       'LABEL'         => $accession,
       'SEQUENCE'      => $extracted_record->{'sequence'},
       'SEQUENCE_TYPE' => 'peptide',
       'SOURCE_ID'     => $self->_get_source_id(),
       'SPECIES_ID'    => $self->{'species_id'},
       'STATUS'        => 'experimental',     # FIXME: seems at least some TrEMBL entries should be 'predicted' instead
       'SYNONYMS'      => \@synonyms,
       '_multiplicity' => $xref_multiplicity, # hint for Loader
     };
  # FIXME: still to be defined: DIRECT_XREFS, DEPENDENT_XREFS
  # Remember not to add UniProt Gene Name dependent xrefs for proteins
  # derived from Ensembl.

  return $xref_graph_node;
}


# FIXME: description
sub get_source_id_map {
  my ( $self ) = @_;

  # Just in case, even though we presently call _load_maps() in the
  # constructor so it shouldn't be possible to call this method before
  # maps have been loaded
  if ( ! exists $self->{'maps'}->{'named_source_ids'} ) {
    croak 'Source-ID map is missing';
  }

  return $self->{'maps'}->{'named_source_ids'};
}


# FIXME: at present these are ALL mock entries
sub _load_maps {
  my ( $self ) = @_;

  # FIXME: $self->baseParserInstance->get_taxonomy_from_species_id( $species_id, $dbh );
  my $mock_taxonomy_ids_for_species
    = {
       9606 => 1, ## no critic(ProhibitMagicNumbers)
     };
  $self->{'maps'}->{'taxonomy_ids_for_species'}
    = $mock_taxonomy_ids_for_species;

  # FIXME: $self->baseParserInstance->get_source_id_for_source_name( $name, $priority, $dbh );
  my $mock_source_id_map
    = {
       'Uniprot/SWISSPROT'
       => {
           'direct'          => 138,
           'sequence_mapped' => 139,
         },
       'Uniprot/SPTREMBL'
       => {
           'direct'            => 134,
           "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD" => 136,
           'sequence_mapped'   => 135,
         },
     };
  $self->{'maps'}->{'named_source_ids'}
    = $mock_source_id_map;

  # FIXME: either abort on any loading failure or implement lazy loading

  return;
}


# Returns true if the current record describes an entry derived from
# Ensembl, false otherwise.
sub _entry_is_from_ensembl {
  my ( $self ) = @_;

  # As of end of October 2018, the old parser's way of identifying
  # proteins derived from Ensembl by searching for a comment topic
  # "CAUTION" stating "The sequence shown here is derived from an
  # Ensembl" no longer works because there seem to be no entries
  # containing this string any more. Therefore, for the time being
  # this test always returns false.

  return 0;
}


# Returns true if the current record describes an entry tagged as
# unreviewed, false otherwise.
sub _entry_is_unreviewed {
  my ( $self ) = @_;

  # This is the way the old UniProtParser identified unreviewed
  # entries. FIXME: is this still a thing? As of October 2018 there
  # are NO such entries in either the first ~1000 lines of the TrEMBL
  # file or anywhere in the SwissProt one.
  my $accession_numbers = $self->{'extracted_record'}->{'accession_numbers'};
  if ( lc( $accession_numbers->[0] ) eq 'unreviewed' ) {
    return 1;
  }

  return 0;
}


# Translate quality of the extracted entry into the matching Ensembl
# source_id.
sub _get_source_id {
  my ( $self ) = @_;

  my $source_id_map = $self->{'maps'}->{'named_source_ids'};

  my $entry_quality = $self->{'extracted_record'}->{'quality'};
  my $criteria = $source_selection_criteria_for_status{ $entry_quality->{'status'} };
  my $priority_mapper = $criteria->[1];

  my ( $source_name, $priority )
    = (
       $criteria->[0],
       $priority_mapper->( $entry_quality->{'evidence_level'} ),
     );

  return $source_id_map->{$source_name}->{$priority};
}


# Translate extracted taxon codes into Ensembl taxonomy IDs, then
# return the number of taxons matching the species ID under
# consideration.
sub _recognised_taxon_ids {
  my ( $self ) = @_;

  my $taxon_codes = $self->{'extracted_record'}->{'taxon_codes'};
  my $tid4s_map = $self->{'maps'}->{'taxonomy_ids_for_species'};

  my @taxonomy_ids;
  foreach my $taxon ( @{ $taxon_codes } ) {
    my $code_mapper
      = $taxonomy_ids_from_taxdb_codes{ $taxon->{'db_qualifier'} };
    push @taxonomy_ids, $code_mapper->( $taxon->{'taxon_code'} );
  }

  my $recognised_taxonomy_ids = 0;
  foreach my $taxonomy_id ( @taxonomy_ids ) {
    $recognised_taxonomy_ids += ( exists $tid4s_map->{$taxonomy_id} );
  }

  return $recognised_taxonomy_ids;
}


1;
