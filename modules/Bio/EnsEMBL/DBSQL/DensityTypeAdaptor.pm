#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::DensityTypeAdaptor
#
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::DensityTypeAdaptor

=head1 SYNOPSIS

my $density_type_adaptor = $database_adaptor->get_DensityTypeAdaptor();
@density_types = @{$density_type_adaptor->fetch_all()};

my $dt = $density_type_adaptor->fetch_by_dbID(12);


=head1 DESCRIPTION

DensityTypeAdaptor - Performs database interaction for DensityType objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::DensityTypeAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [1]    : see superclass (Bio::EnsEMBL::DBSQL::BaseAdaptor) arguments
  Example    : #use this instead of the constructor directly:
               my $dta = $db_adaptor->get_DensityTypeAdaptor();
  Description: Constructor. Creates a new DensityTypeAdaptor
  Returntype : Bio::EnsEMBL::DBSQL::DensityTypeAdaptor
  Exceptions : none
  Caller     : DBAdaptor

=cut

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->{'dbID_cache'} = {};

  return $self;
}



=head2 fetch_all

  Arg [1]    : none
  Example    : my @density_types = @{$density_type_adaptor->fetch_all};
  Description: Retrieves every density type in the database
  Returntype : reference to list of Bio::EnsEMBL::DensityType objects
  Exceptions : none
  Caller     : general, new

=cut

sub fetch_all {
  my $self = shift;

  my @out;

  my $sth = $self->prepare("SELECT density_type_id, analysis_id, block_size,".
                           "       value_type " .
                           "FROM density_type");

  $sth->execute();

  my($dbID, $analysis_id, $blk_size, $vtype);
  $sth->bind_columns(\$dbID, \$analysis_id, \$blk_size, \$vtype);

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();

  while($sth->fetch()) {
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);

    my $dt = Bio::EnsEMBL::DensityType->new(-ADAPTOR => $self,
                                            -DBID    => $dbID,
                                            -ANALYSIS => $analysis,
                                            -BLOCK_SIZE => $blk_size,
                                            -VALUE_TYPE => $vtype);

    $self->{'dbID_cache'}->{$dbID} = $dt;

    push @out, $dt;
  }

  return \@out;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : my $dt = $density_type_adaptor->fetch_by_dbID($dbID);
  Description: Retrieves a density type object via its internal identifier
  Returntype : Bio::EnsEMBL::DensityType
  Exceptions : throw if dbID argument not defined
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!defined($dbID)) {
    throw("dbID argument must be defined");
  }

  if($self->{'dbID_cache'}->{$dbID}) {
    return $self->{'dbID_cache'}->{$dbID};
  }

  # go back to database and refill caches
  $self->fetch_all();

  return $self->{'dbID_cache'}->{$dbID};
}


=head2 fetch_all_by_logic_name

  Arg [1]    : string $logic_name
  Example    : my @dts = @{$dtype_adaptor->fetch_all('repeat_coverage')};
  Description: Retrieves all density types with a given logic name
  Returntype : reference to list of Bio::EnsEMBL::DensityTypes
  Exceptions : thrown if logic_name argument is not provided
  Caller     : general

=cut

sub fetch_all_by_logic_name {
  my $self = shift;
  my $logic_name = shift;

  if(!defined($logic_name)) {
    throw("logic_name argument is required.");
  }

  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();

  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  return [] if(!$analysis);

  my $sth = $self->prepare("SELECT density_type_id, block_size,".
                           "       value_type " .
                           "FROM density_type " .
                           "WHERE analysis_id = ?");
  $sth->execute($analysis->dbID());

  my($dbID, $blk_size, $vtype);
  $sth->bind_columns(\$dbID, \$blk_size, \$vtype);

  my @out;

  while($sth->fetch()) {
    my $dt = Bio::EnsEMBL::DensityType->new(-ADAPTOR => $self,
                                            -DBID    => $dbID,
                                            -ANALYSIS => $analysis,
                                            -BLOCK_SIZE => $blk_size,
                                            -VALUE_TYPE => $vtype);

    $self->{'dbID_cache'}->{$dbID} = $dt;

    push @out, $dt;
  }

  return \@out;
}


#
# garbage collection method, automatically called when DBAdaptor is cleaned up
#
sub deleteObj {
  my $self = shift;

  delete $self->{'dbID_cache'};

  $self->SUPER::deleteObj;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DensityType @dt
               the density types to store in the database
  Example    : $density_type->store(1234, @density_types);
  Description: Stores a list of density type objects in the database
  Returntype : none
  Exceptions : thrown if @dt is not defined
               or if any elements of @dt are not Bio::EnsEMBL::DensityType 
  Caller     : general

=cut

sub store {
  my ($self,@dt) = @_;

  if( scalar(@dt) == 0 ) {
    throw("Must call store with list of Density Types");
  }

  my $sth = $self->prepare
    ("INSERT INTO density_type (analysis_id,".
                                  "block_size, value_type) ". 
    "VALUES (?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $dt ( @dt ) {
    if( !ref $dt || !$dt->isa("Bio::EnsEMBL::DensityType") ) {
      throw("Density Type must be an Ensembl DensityType, " .
            "not a [".ref($dt)."]");
    }

    if($dt->is_stored($db)) {
      warning("DensityType [".$dt->density_type_id."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($dt->analysis())) {
      throw("An analysis must be attached to the density type to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$dt->analysis->is_stored($db)) {
      $analysis_adaptor->store($dt->analysis());
    }
	
    $sth->execute($dt->analysis->dbID(), $dt->block_size(), $dt->value_type());

    my $dbID = $sth->{'mysql_insertid'};

    # next two lines are to set the density type as stored
    $dt->dbID($dbID);
    $dt->adaptor($self);

    $self->{'dbID_cache'}->{$dbID} = $dt;
  }
}

1;
