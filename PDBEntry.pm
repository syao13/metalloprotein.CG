## PDBEntry.pm
#
package PDBEntry;
use strict;
use Atom;
use Sequence;
use Residue;

our @defaultDataMembers = ("singlePdbFile" => "" );

## usage PDBEntry->new(singlePdbFile => "file_path")
sub new
  {
  my $class = shift @_;

#  my $self = { @defaultDataMembers, @_, "zn" => $zn, "atoms" => 0 };
  my $self = { @defaultDataMembers, "atoms" => [], "residues" => [], "_residueIDmap" => {}, "sequences" => {}, @_ };
  bless $self, ref $class || $class;
  
  $self->read() if ($self->{singlePdbFile} ne "");
  return $self;  
  }

##
sub read
  {
  my $self = shift @_;
  my $singlePdbFile = (@_) ? shift @_ : $self->{singlePdbFile};
  $self->{singlePdbFile} = $singlePdbFile;
  my $atoms = (@_) ? shift @_ : [];
  $self->{atoms} = $atoms;
  $self->{residues} = [];
  $self->{_residueIDmap} = {};

  if ($singlePdbFile =~ /\.gz$/)
    { open (PDBFILE, "/bin/zcat $singlePdbFile|") || die "Error in opening $singlePdbFile."; }
  else
    { open (PDBFILE, "<$singlePdbFile") || die "Error in opening $singlePdbFile."; }

  my $modelCount = 0;
  my $method;
  my $PDBid;
  my $date;
  my $chainSeqs = {};
  my $resolution = -1;
  my $rValue = -1;
  my $rFree = -1;

  foreach my $record (<PDBFILE>)
    {
    chomp $record;
    if ($record =~ /^ATOM/ || $record =~ /^HETATM/)
      {
      my @keyValues = ('record' => $record,
                        'recordType' => substr($record, 0, 6),
                        'serial' => substr($record, 6, 5), 
                        'atomName'=> substr($record, 12, 4),
                        'residueName'=> substr($record, 17, 3), 
                        'chainID' => substr($record, 21, 1),
                        'residueNumber' => substr($record, 22, 4), 
                        'x' => substr($record, 30, 8),
                        'y' => substr($record, 38, 8), 
                        'z' => substr($record, 46, 8),
			'occupancy' => substr($record, 54, 6),
			'bFactor' => substr($record, 60, 6),  
                        'element' => substr($record, 76, 2),
                        'file' => $singlePdbFile, 
                        'PDBid' => $PDBid, 
                        'method' => $method,
			'date' => $date,
			'resolution' => $resolution,
			'rValue' => $rValue,
			'rFree' => $rFree
                       ) ;
      foreach my $value (@keyValues)
        { 
        $value =~ s/^\s+//; ## Remove leading white spaces.
        $value =~ s/\s+$//; ## Remove tailing white spaces.
        }
      
      push @$atoms, Atom->new(@keyValues);
      }
    elsif ($record =~ /^HEADER/)
      { 
      $PDBid = substr($record, 62, 4); 
      $date = substr($record, 57, 2);
      }
    elsif ($record =~ /^EXPDTA/)
      { 
      $method = substr($record, 6, 30);
      $method =~ s/^\s+//; ## Remove leading white spaces.
      $method =~ s/\s+$//; ## Remove tailing white spaces.
      $method =~ s/\s/_/g;
      }
    elsif ($record =~ /^REMARK   2 RESOLUTION.(.+)ANGSTROMS/)
      { $resolution = $1;}
    elsif ($record =~ /^REMARK   3   R VALUE            \(WORKING SET\) : (.+)$/)
      { $rValue = $1; }
    elsif ($record =~ /^REMARK   3   FREE R VALUE                     : (.+)$/)
      {$rFree = $1; }
    elsif ($record =~ /^MODEL/)
      {
      $modelCount++;
      last if ($modelCount > 1); ## Use the first model if there are more the one model structure
      }
    elsif ($record =~ /^SEQRES/)
      {
      my $chain = substr($record, 11, 1);
      my $chainSeqOne = substr($record, 19, 51);

      push (@{$$chainSeqs{$chain}}, $chainSeqOne);
      }
    
    }
  close(PDBFILE);

  my $residueMap = $self->{_residueIDmap};
  foreach my $atom (@$atoms)
    {
    my $resID = $atom->resID();
    if (! exists $$residueMap{$resID})
      { $$residueMap{$resID} = Residue->new($atom); }
    else
      { $$residueMap{$resID}->add($atom); }
    }  

  @{$self->{residues}} = sort { $a->cmp($b); } (values %$residueMap);

  my $seqs = [];      
  foreach my $chain (keys %$chainSeqs)
    {
    push (@$seqs, Sequence->createSEQRESseq($chain, $$chainSeqs{$chain}, $atoms)) ;
    }
  $self->{sequences}->{seqres} = $seqs; #SEQRES sequences
  $self->{sequences}->{atom} = Sequence->createATOMseqs($atoms); #ATOM sequences
  $self->{sequences}->{number} = Sequence->createATOMnum($atoms); #ATOM numbers
 
  return $self;
  }


sub residue
  {
  my $self = shift @_;
  my $resID = shift @_;
  
  return (exists $self->{_residueIDmap}{$resID}) ? $self->{_residueIDmap}{$resID} : 0;
  }


1;
