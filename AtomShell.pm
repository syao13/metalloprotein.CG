## AtomShell.pm 
#
package AtomShell;
use strict;
use Atom;

## Constructor: creatShells; creat; new
## Methods: metalID, anglesBetweenShell, distanceBetweenShell

our @defaultDataMembers = (   
                                "center" => 0,
                                "shell" => 0,
                                "secondShell" => 0,
				"seqsOfPDB" => 0
);

our $secondShellDist = 8.0;
our %atomRadius = (
	"ZN" =>	1.35,
	"YB" => 1.75,
	"Y"  =>	1.80,
	"W"  =>	1.35,
	"V"  =>	1.35,
	"TL" =>	1.90,
	"TB" =>	1.75,
	"SR" =>	2.00,
	"SM" =>	1.85,
	"SB" =>	1.45,
	"RU" =>	1.30,
	"RH" =>	1.35,
	"RB" =>	2.35,
	"PT" =>	1.35,
	"PR" =>	1.85,
	"PD" =>	1.40,
	"PB" =>	1.80,
	"OS" =>	1.30,
	"NI" =>	1.35,
	"NA" =>	1.80,
	"MO" =>	1.45,
	"MN" =>	1.40,
	"MG" =>	1.50,
	"LU" =>	1.75,
	"LI" =>	1.45,
	"LA" =>	1.95,
	"K"  => 2.20,
	"IR" =>	1.35,
	"IN" =>	1.55,
	"HO" =>	1.75,
	"HG" =>	1.50,
	"GD" =>	1.80,
	"GA" =>	1.30,
	"FE" =>	1.40,
	"EU" =>	1.85,
	"ER" =>	1.75,
	"DY" =>	1.75,
	"CU" =>	1.35,
	"CS" =>	2.60,
	"CR" =>	1.40,
	"CO" =>	1.35,
	"CE" => 1.85,
	"CD" =>	1.55,
	"CA" =>	1.80,
	"BI" =>	1.60,
	"BA" =>	2.15,
	"AU" => 1.35,
	"AL" => 1.25,
	"AG" => 1.60);


## Creat a list of AtomShell objs of given element from a list of Atom objects 
## Optional parameters: min and max distances, with default as 1.3 to 3.2
sub createShells
  {
  my $class = shift @_;  $class = ref $class || $class;
  my $element = shift @_;
  my $atoms = shift @_;
  my $minDist = (@_) ? shift @_ : 1.3;

  my $maxDist = (@_) ? shift @_ : (($atomRadius{$element} > $atomRadius{"ZN"})? ($atomRadius{$element} - $atomRadius{"ZN"} + 3.2): 3.2 );
  my $ligElements = (@_) ? shift @_ : ();

  my @centers = grep { $_->{element} eq $element && substr($_->{chainID}, 0,1) ne "#"; } (@$atoms); 

  ## Rule out clusters
  my %cluster;
  foreach my $i (0..(@centers-1))
    {
    foreach my $j (($i+1)..(@centers))
      {
      if ($centers[$i]->distance($centers[$j]) < 3)
        {
        $cluster{$centers[$i]} = 1;
        $cluster{$centers[$j]} = 1;
	}
      }
    }
  @centers = (grep {! exists $cluster{$_}; } (@centers));
  #print scalar @centers, "\n";

  return [ map { $class->create($_,$atoms,$minDist,$maxDist,$ligElements); } (@centers) ];
  }

## create an AtomShell obj from a center atom and list of Atom objs.
sub create
  {
  my $class = shift @_;  $class = ref $class || $class;
  my $center = shift @_;
  my $atoms = shift @_;
  my $minDist = shift @_;
  my $maxDist = shift @_;
  my $ligElements = shift @_;

  ## if elements are specified, use them; otherwise using everything other than H.
  my (@ligElements, @tempShell);
  if ($ligElements)
    {
    @ligElements = split(/(?=[A-Z]+[^A-Z]?)/, $ligElements);
    @tempShell = grep {my $distance = $center->distance($_); my $atom = $_; ($distance >= $minDist && $distance <= $maxDist && grep {$atom->{element} eq uc($_)} (@ligElements) );} (@$atoms);
#print $center->{PDBid}, ".", $center->{chainID}, ".", $center->{residueNumber}, "\n" if (grep {$center->distance($_) < $minDist && $center->distance($_) > 0} (@$atoms));
    }
  else 
    { @tempShell = grep {my $distance = $center->distance($_); ($distance >= $minDist && $distance <= $maxDist && $_->{element} ne "H" && ! $atomRadius{$_->{element}} );} (@$atoms); }
#print $center->coordinates(), ", center\n";
#map {print $_->coordinates(), ", ligands\n";} (grep {$_->{residueNumber} == 126} (@$atoms));

#print $center->resID(), "---before\n";
#map {print $_->resID(), "\n";} (@tempShell);

  my @remove;
  foreach my $atom (@tempShell)
    {
    if (substr ($atom->{chainID}, 0, 1) eq "#") 
      {
      foreach my $rest (grep {$_ ne $atom} (@tempShell))
        {
        if ($rest->distance($atom) < 0.00001)
	  { push @remove, $atom; }
        }
      }
    }
  @tempShell = grep {my $atom = $_; grep {$_ ne $atom} (@remove) ;} (@tempShell) if (@remove);

print $center->{PDBid}, ".", $center->{chainID}, ".", $center->{residueNumber}, "\n" if (@remove);
#print $center->resID(), "--- after\n";
#map {print $_->resID(), "\n";} (@tempShell);

  return $class->new("center" => $center, 
		     "shell" => [_remove2ndShell($center, @tempShell)], 
		     "secondShell" => [ grep {my $distance = $center->distance($_); ($distance >= $minDist && $distance <= $secondShellDist && $_->{element} ne "C" && $_->{element} ne "H" && ! $atomRadius{$_->{element}} );} (@$atoms) ] ); 
  }


sub _remove2ndShell
  {
  my $center = shift @_;

  my %remove;
  foreach my $ligand (@_)
    {
    if (grep {$center->distance($ligand) > $center->distance($_) * 1.5 && $center->distance($ligand) > $ligand->distance($_) * 1.5 ;} (@_))
      { 
      $remove{$ligand} = 1;

      #my $pdbid = $center->{PDBid};
      #my $chainid = $center->{chainID}; 
      #my $serial = $center->{residueNumber};
      #my $residue = $ligand->resID();
      #my $ele = $ligand->{element};

      # if elements are specified, then there normally shouldn't be any C; and if not specified, then we want to include C to estimate the density issue.
      #if ($ele ne "C")
	#{
	#print  "$pdbid.$chainid.$serial:$residue.$ele\n";
        #print join (", ", map {$_->resID(), $_->{element}, $center->distance($_)} (@_)), "\n";
	#}
      }
    }
  
  return grep { ! $remove{$_} ; } (@_);  
  }

## A standard way to creat a now obj, not useful for this AtomShell obj though
sub new
  {
  my $class = shift @_;
  my $self = { @defaultDataMembers, @_ };
  bless $self, ref $class || $class;

  return $self;  
  }

## angles_between_atoms
##   calculates every atom to atom angle in a set of atoms
##   &angles_between_atoms ( ref to array of atom information hashes, zinc information hash)
sub anglesBetweenShell
  {
  my $self = shift @_;

  my @angles;
  foreach my $a (0..($#{$self->{shell}} - 1))
    { push @angles, map { $self->{center}->angle($self->{shell}[$a], $self->{shell}[$_]); } (($a+1)..($#{$self->{shell}})); }

  return @angles;
  }

## distance_between_atoms
##    calculates pairwise distances between a list of atoms
## Parameters:
##    $atoms - reference to array of atom info hash references
sub distanceBetweenShell
  {
  my $self = shift @_;
  my $atoms = shift @_;
  my @distances;

  foreach my $a (0..($#{$self->{shell}} - 1))
    { push @distances, map {$self->{shell}[$a]->distance($self->{shell}[$_]);} ( ($a+1)..($#{$self->{shell}}) ) ; }

  return @distances;
  }


## Returns a zinc id that could be used to identify each zinc site specifically
## 	The format is pdbid.chinid.serial#, for example, 1BY4.C.2330
sub metalID
  {
  my $self = shift @_;

  my $metalAtom = $self->{center};

  my $pdbid = $metalAtom->{PDBid};
  my $chainid = $metalAtom->{chainID};
  my $serial = $metalAtom->{residueNumber};

  return  "$pdbid.$chainid.$serial";
  }



1;
