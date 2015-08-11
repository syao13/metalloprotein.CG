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
  my $maxDist = ($atomRadius{$element} > $atomRadius{"ZN"})? ($atomRadius{$element} - $atomRadius{"ZN"} + 3.2): 3.2;

  my @centers = grep { $_->{element} eq $element; } (@$atoms); 

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

  return [ map { $class->create($_,$atoms,$minDist,$maxDist); } (@centers) ];
  }

## create an AtomShell obj from a center atom and list of Atom objs.
sub create
  {
  my $class = shift @_;  $class = ref $class || $class;
  my $center = shift @_;
  my $atoms = shift @_;
  my $minDist = shift @_;
  my $maxDist = shift @_;

  #print $center->{element}, "flag\n";
  return $class->new("center" => $center, 
		     "shell" => [ grep {my $distance = $center->distance($_); ($distance >= $minDist && $distance <= $maxDist && $_->{element} ne "C" && $_->{element} ne "H" && ! $atomRadius{$_->{element}} );} (@$atoms) ],
		     "secondShell" => [ grep {my $distance = $center->distance($_); ($distance >= $minDist && $distance <= $secondShellDist && $_->{element} ne "C" && $_->{element} ne "H" && ! $atomRadius{$_->{element}} );} (@$atoms) ] ); 
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
