#!/usr/local/ensembl/bin/perl

=head1 NAME

check_mapping.pl - script to check whole genome alignment between two
assemblies.

=head1 SYNOPSIS

check_mapping.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY

  --altdbname=NAME                    alternative database NAME
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

  --althost=HOST                      alternative databases host HOST
  --altport=PORT                      alternative database port PORT
  --altuser=USER                      alternative database username USER
  --altpass=PASS                      alternative database password PASS

  --chromosomes, --chr=LIST           only process LIST toplevel seq_regions

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script checks if the whole genome alignment between two assemblies is
correct. It does so by comparing the sequence in the reference database with
the sequence of the projected fragments in the alternative database.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altdbname=s',
    'altassembly=s',
    'althost=s',
    'altport=n',
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'althost',
    'altport',
    'chromosomes',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
    'assembly',
    'altdbname',
    'altassembly'
);

# first set connection parameters for alternative db if not different from
# reference db
map { $support->param("alt$_", $support->param($_)) unless ($support->param("alt$_")) } qw(host port user);

# reference database
my $R_dba = $support->get_database('ensembl');
my $R_sa = $R_dba->get_SliceAdaptor;

# database containing the alternative assembly
my $A_dba = $support->get_database('core', 'alt');
my $A_sa = $A_dba->get_SliceAdaptor;

$support->log("Looping over toplevel seq_regions...\n\n");

foreach my $chr ($support->sort_chromosomes) {
  $support->log_stamped("Toplevel seq_region $chr...\n", 1);

  my $R_slice = $R_sa->fetch_by_region('toplevel', $chr);
  my $A_slice = $A_sa->fetch_by_region('toplevel', $chr);

  unless ($A_slice) {
    $support->log("Not found in alternative db. Skipping.\n", 2);
    next;
  }
  
  my $cs_name = $A_slice->coord_system_name;

  # compare reference and alternative sequence
  my @segments = @{ $R_slice->project($cs_name, $support->param('altassembly')) };
  
  my $i;
  my $k;

  foreach my $seg (@segments) {
    # reference sequence
    my $R_sub_slice = $R_slice->sub_Slice($seg->from_start, $seg->from_end);
    my $R_seq = $R_sub_slice->seq;
    
    # alternative sequence
    my $A_proj_slice = $seg->to_Slice;
    
    # ignore PAR region (i.e. we project to the symlinked seq_region)
    next if ($A_proj_slice->seq_region_name ne $chr);
    
    my $A_sub_slice = $A_slice->sub_Slice($A_proj_slice->start, $A_proj_slice->end, $A_proj_slice->strand);
    my $A_seq = $A_sub_slice->seq;

    # compare
    if ($R_seq eq $A_seq) {
      # sequences are identical -> ok
      $support->log_verbose("Sequence match at ".$R_sub_slice->name."\n", 2);

    } else {
      # not identical -> something is wrong
      $support->log("Sequence mismatch at ".$R_sub_slice->name."\n", 2);

      my $R_sub_seq;
      my $A_sub_seq;
      
      if ($R_sub_slice->length > 20) {
        $R_sub_seq = substr($R_seq, 0, 10)."...".substr($R_seq, -10, 10);
        $A_sub_seq = substr($A_seq, 0, 10)."...".substr($A_seq, -10, 10);
      } else {
        $R_sub_seq = substr($R_seq, 0, 20);
        $A_sub_seq = substr($A_seq, 0, 20);
      }
      
      $support->log("Ref: $R_sub_seq\n", 3);
      $support->log("Alt: $A_sub_seq\n\n", 3);

      $i++;
    }

    $k++;
  }

  if ($i) {
    $support->log("Total: $i (of $k) alignments contain sequence mismatches.\n", 2);
  } else {
    $support->log("All $k alignments ok.\n", 2);
  }

  $support->log_stamped("Done.\n\n", 1);
}

# finish logfile
$support->finish_log;

