## PDBEntry.pm
#
package PDBEntry;
use strict;
use Atom;
use Sequence;
use Residue;
use Math::MatrixReal;

our @defaultDataMembers = ("singlePdbFile" => "" , "metal" => "");
our %metalFull = ("ZN" => "ZINC", "FE" => "IRON", "CA" => "CALCIUM", "NA" => "SODIUM", "MG" => "MAGNESIUM");

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
  my ($method, $PDBid, $date);
  my $chainSeqs = {};
  my $resolution = -1;
  my $rValue = -1;
  my $rFree = -1;
  my $solvent = "False";
  my $crystalMats = [];
  my $bioMats = [];
  my $ortha = [];
  my $orthb = [];
  my $orthc = [];
  my $biomolecules = [];
  my $authorU = 0;

  while ( my $record = <PDBFILE>)
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
			'rFree' => $rFree,
			'solvent' => $solvent
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
    elsif ($record =~ /^REMARK 280.* ($self->{metal}|$metalFull{$self->{metal}})/)
      {$solvent = "True";}
    elsif ($record =~ /^REMARK 290   SMTRY/)
      {
      my @mat = split (/\s+/, $record);
      my $smtryNum = substr ($mat[2], -1); 
      $$crystalMats[$mat[3]-1][$smtryNum] = [0, $mat[4], $mat[5], $mat[6]]; ## add 0 in front so that the index will match
      $$crystalMats[$mat[3]-1][4][$smtryNum] = $mat[7]; 
      }
    elsif ($record =~ /^REMARK 350 BIOMOLECULE:/)
      {
      ## Use author-determined instead of sofeware-determined biological symmetry matrix if there is any. Otherwise use the sofeware-determined.
      my $authorL;
      my $bioMats = [];
      my $matChains = [];
      while (my $nextline = <PDBFILE>) 
	{
        if ($nextline =~ /^REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT/)
	  {
	  $authorL = 1;
	  $authorU = 1; 
	  }
	elsif ($nextline =~ /^REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE/)
	  {
	  next if $authorL == 1;  
	  last if $authorU == 1; 
	  }
	elsif ($nextline =~ /^REMARK 350 APPLY THE FOLLOWING TO CHAINS: (.+)$/)
	  {
	  my $var = $1;
	  $var =~ s/\s+$//; 
          $var =~ s/,$// if substr($var, -1) eq ",";
	  @$matChains = split (', ', $var); 
	  }
        elsif ($nextline =~ /^REMARK 350                    AND CHAINS: (.+)$/)
	  {
	  my $var = $1;
          $var =~ s/\s+$//;
          $var =~ s/,$// if substr($var, -1) eq ","; 
	  push @$matChains, (split (', ', $var));
	  } 
	elsif ($nextline =~ /^REMARK ...                                     /)
	  { last; }
	elsif ($nextline =~ /^REMARK 350   BIOMT/)
	  {
	  my @mat = split (/\s+/, $nextline);
          my $smtryNum = substr ($mat[2], -1);
          $$bioMats[$mat[3]-1][$smtryNum] = [0, $mat[4], $mat[5], $mat[6]]; ## add 0 in front so that the index will match
          $$bioMats[$mat[3]-1][4][$smtryNum] = $mat[7];
	  }
	}
      push @$biomolecules, {"matrices" => $bioMats, "chains" => $matChains} ;
      }
    elsif ($record =~ /^SEQRES/)
      {
      my $chain = substr($record, 11, 1);
      my $chainSeqOne = substr($record, 19, 51);

      push (@{$$chainSeqs{$chain}}, $chainSeqOne);
      }
    elsif ($record =~ /^CRYST1/)
      { my @cryst = split (/\s+/, $record); }
    elsif ($record =~ /^SCALE1/)
      { $ortha = [(split (/\s+/, $record))[1..3]]; }
    elsif ($record =~ /^SCALE2/)
      { $orthb = [(split (/\s+/, $record))[1..3]]; }
    elsif ($record =~ /^SCALE3/)
      { $orthc = [(split (/\s+/, $record))[1..3]]; }
    }
  close(PDBFILE);

  ## comment out residue calc for now, 2016. 3.13, as no longer calc chi dihedreal angles.
  #my $residueMap = $self->{_residueIDmap};
  #foreach my $atom (@$atoms)
  #  {
  #  my $resID = $atom->resID();
  #  if (! exists $$residueMap{$resID})
  #    { $$residueMap{$resID} = Residue->new($atom); }
  #  else
  #    { $$residueMap{$resID}->add($atom); }
  #  }  

  #@{$self->{residues}} = sort { $a->cmp($b); } (values %$residueMap);

  my $seqs = [];      
  foreach my $chain (keys %$chainSeqs)
    {
    push (@$seqs, Sequence->createSEQRESseq($chain, $$chainSeqs{$chain}, $atoms)) ;
    }
  $self->{sequences}->{seqres} = $seqs; #SEQRES sequences
  $self->{sequences}->{atom} = Sequence->createATOMseqs($atoms); #ATOM sequences
  $self->{sequences}->{number} = Sequence->createATOMnum($atoms); #ATOM numbers
 
  goto THEEND if $method ne "X-RAY_DIFFRACTION";

  ## biological symmetry
  my @bioAtoms;
  foreach my $bio (@$biomolecules)
    {
    next if scalar @{$$bio{matrices}} == 1;

    foreach my $num (1..(@{$$bio{matrices}}-1))
      {
      my @atomsChain = grep {my $atom = $_; grep {$_ eq $atom->{chainID};} (@{$$bio{chains}}); } (@$atoms);
      my @atomsNew = map { $_->transformBio($$bio{"matrices"}[$num]); } (@atomsChain);
      push @bioAtoms, @atomsNew;      
      }
    }

  ## symmetry related atoms
  ## x' = Rx + T 
  my $xmin = (sort {$a <=> $b} (map {$_->{x}} (@$atoms)))[0];
  my $xmax = (sort {$b <=> $a} (map {$_->{x}} (@$atoms)))[0];
  my $ymin = (sort {$a <=> $b} (map {$_->{y}} (@$atoms)))[0];
  my $ymax = (sort {$b <=> $a} (map {$_->{y}} (@$atoms)))[0];
  my $zmin = (sort {$a <=> $b} (map {$_->{z}} (@$atoms)))[0];
  my $zmax = (sort {$b <=> $a} (map {$_->{z}} (@$atoms)))[0];

  ## Get orthogonalization matrix O from deororthogonalization matrix O' (given by SCALE record)
  ## The neighbering cells can be calculated using formula, X' = O(O'(RX + T) + T') = OO'(RX+T) + OT' = RX+T + O[-1/0/1,-1/0/1,-1/0/1] 
  my $inverseO = Math::MatrixReal->new_from_rows([$ortha, $orthb, $orthc]);
  my $bigO = $inverseO->inverse();

  my @symAtoms;
  foreach my $i (-1, 0, 1)
    {
    foreach my $j (-1, 0, 1)
      {
      foreach my $k (-1, 0, 1)
	{
	foreach my $t (0..(@$crystalMats-1))	
	  {
	  next if ($i == 0 && $j == 0 && $k == 0 && $t == 0);
	  my $mat = $$crystalMats[$t];

	  my $neighborX = $bigO->element(1,1) * $i + $bigO->element(1,2) * $j + $bigO->element(1,3) * $k;
          my $neighborY = $bigO->element(2,1) * $i + $bigO->element(2,2) * $j + $bigO->element(2,3) * $k;
          my $neighborZ = $bigO->element(3,1) * $i + $bigO->element(3,2) * $j + $bigO->element(3,3) * $k;	

	  my @xs = ($$mat[1][1] * $xmin + $$mat[1][2] * $ymin + $$mat[1][3] * $zmin + $$mat[4][1] + $neighborX, $$mat[1][1] * $xmax + $$mat[1][2] * $ymax + $$mat[1][3] * $zmax + $$mat[4][1] + $neighborX) ;
          my @ys = ($$mat[2][1] * $xmin + $$mat[2][2] * $ymin + $$mat[2][3] * $zmin + $$mat[4][2] + $neighborX, $$mat[2][1] * $xmax + $$mat[2][2] * $ymax + $$mat[2][3] * $zmax + $$mat[4][2] + $neighborX) ;
          my @zs = ($$mat[3][1] * $xmin + $$mat[3][2] * $ymin + $$mat[3][3] * $zmin + $$mat[4][3] + $neighborX, $$mat[3][1] * $xmax + $$mat[3][2] * $ymax + $$mat[3][3] * $zmax + $$mat[4][3] + $neighborX) ;

	  @xs = sort {$a <=> $b} (@xs);
          @ys = sort {$a <=> $b} (@ys);
          @zs = sort {$a <=> $b} (@zs);

  	  if (($xs[0] < $xmax + 3.5) && ($xs[1] > $xmin - 3.5) && ($ys[0] < $ymax + 3.5) && ($ys[1] > $ymin - 3.5)  && ($zs[0] < $zmax + 3.5) && ($zs[1] > $zmin - 3.5) )
	    {
	    my @atomsNew = map { $_->transform($mat, $neighborX, $neighborY, $neighborZ); } (@$atoms);
            push @symAtoms, @atomsNew;
	    }

#map { print join (", ", @$_), "\n" if $_;} (@$mat);
#my $atom = (grep {$_->{residueNumber} == 126 && $_->{atomName} eq "OE1"} (@$atoms))[0];
#my $newAtom = $atom->transform($mat, $i * $crystA, $j * $crystB, $k * $crystC);
#print $atom->coordinates(), "\n";
#print $newAtom->coordinates(), "\n";
#print "originals, $xmin, $xmax, $ymin, $ymax, $zmin, $zmax\n";
#print "new, $xminNew, $xmaxNew, $yminNew, $ymaxNew, $zminNew, $zmaxNew\n";
	  }
	}
      }
    }
  push @$atoms, @symAtoms;

  THEEND:
  return $self;
  }


sub residue
  {
  my $self = shift @_;
  my $resID = shift @_;
  
  return (exists $self->{_residueIDmap}{$resID}) ? $self->{_residueIDmap}{$resID} : 0;
  }


1;
