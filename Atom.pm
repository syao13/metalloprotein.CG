## Atom.pm
##
package Atom;
use strict;

## Constructor: new
## Methods: resID, angle, distance, torsionAngle

our @defaultAtom = ( 
                 "atomName" => "",
                 "residueType" => "", 
                 "chainID" => "",
                 "residueNumber" => "", 
                 "x" => 0,
                 "y" => 0, 
                 "z" => 0, 
                 "element" => "",
                 # PDB entry specific members
                 "serial" => 0,
                 "recordType" => "",
                 "record" => "", 
                 "file" => "", 
                 "PDBid" => "", 
                 "method" => ""
                 ); 


sub new
  {
  my $class = shift @_;
  my $self = { @defaultAtom, @_};

  return bless $self, ref $class || $class;
  }


## find the closest aa of given atom from atoms
sub closest
  {
  my $self = shift @_;
  my $atoms = shift @_;

  my $minDist = 99;
  my $closestAtom;
  foreach my $atom (@$atoms)
    {
    next if (! &Sequence::_aaCode($atom->{residueName}));
    next if ($self eq $atom);
    next if ($atom->{element} eq "C" || $atom->{element} eq "H" || $atom->{element} eq "ZN") ;
 
    my $dist = $self->distance($atom);
    if ($dist < $minDist)
      {
      $minDist = $dist;
      $closestAtom = $atom;
      }
    }

  return $closestAtom;
  }

## resID
## 	returns the residue of a given atom
##	Mainly for finding chi-1 torsion angle
sub resID
  { 
  my $self = shift @_;
  return $self->{chainID} . $self->{residueNumber} ; 
  }


##  angle
###   calculates the angle between the given two points in a triangle, given object point
###   uses &distance
###   each point input needs to be a hash that contains x y and z coords
###   &angle( point a information hash, point b information hash, point c information hash )
###         (  atom hash, atom 2 hash, zinc hash )  <---   will calculate atom to atom 2 angle
sub angle
  {
  my $self = shift @_;
  my $a = shift @_;
  my $b = shift @_;

  my $distAB = &distance($a, $b);
  my $distASelf = &distance($a, $self);
  my $distBSelf = &distance($b, $self);
  if ($distAB == 0 || $distASelf == 0 || $distBSelf == 0)
    { return 0; }
  if (abs($distAB - $distASelf - $distBSelf) < 0.00001) 
    { return 180; }
  my $angleAB = 57.2957795*(&_arccos( ($distAB**2 - $distASelf**2 - $distBSelf**2)/-(2*$distASelf*$distBSelf) ));
  return $angleAB;
  }


## distance
###   calculates distance between two atoms
###   &distance (hash reference of first atom, hash reference of second atom)
sub distance
  {
  my $self = shift @_;
  my $b = shift @_;

  my $distance = sqrt( ($self->{x} - $b->{x})**2 + ($self->{y} - $b->{y})**2 + ($self->{z} - $b->{z})**2 );
  return $distance;
  }

## torsionAngle
## Calculate dihedral angle of p1-p2-p3 and p2-p3-p4, consider the corresponding four atom chain as three vectors, b1, b2, b3
## Calculation method is from http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
sub torsionAngle
  {
  my $self = shift @_;
  my $point2 = shift @_;
  my $point3 = shift @_;
  my $point4 = shift @_;

  my $b1 = &_vec2points($self, $point2);
  my $b2 = &_vec2points($point2, $point3);
  my $b3 = &_vec2points($point3, $point4);

  my $cross1 = &_crossProduct($b1, $b2);
  my $cross2 = &_crossProduct($b2, $b3);

  my $n1 = &_norm($cross1);
  my $n2 = &_norm($cross2);
  my $nb2 = &_norm($b2);

  my $m1 = &_crossProduct($n1, $nb2);

  my $x = &_dotProduct($n1, $n2);
  my $y = &_dotProduct($m1, $n2);

  my $angle = -57.2957795*(atan2($y, $x));
  return $angle;
  }

sub perpendicularVec
  {
  my $self = shift @_;
  my $a = shift @_;
  my $b = shift @_;

  my $oa = &_vec2points($self, $a);
  my $ob = &_vec2points($self, $b);

  return &_crossProduct($oa, $ob);
  }


## The rest functions are internal functons for trosion anlge calculation only
sub _vec2points
  {
  my ($point1, $point2) = @_;

  my $vec = [];
  $$vec[0] = $$point2{"x"} - $$point1{"x"};
  $$vec[1] = $$point2{"y"} - $$point1{"y"};
  $$vec[2] = $$point2{"z"} - $$point1{"z"};

  return $vec;
  }


sub _crossProduct
  {
  my ($vec1, $vec2) = @_;

  my $crossVec = [];
  $$crossVec[0] = $$vec1[1] * $$vec2[2] - $$vec1[2] * $$vec2[1];
  $$crossVec[1] = $$vec1[2] * $$vec2[0] - $$vec1[0] * $$vec2[2];
  $$crossVec[2] = $$vec1[0] * $$vec2[1] - $$vec1[1] * $$vec2[0];

  return $crossVec;
  }

sub _norm
  {
  my $vec = shift @_;

  my $length = &_vecLength($vec);

  my $normVec = [];

  $$normVec[0] = $$vec[0] / $length;
  $$normVec[1] = $$vec[1] / $length;
  $$normVec[2] = $$vec[2] / $length;

  return $normVec;
  }

sub _dotProduct
  {
  my ($vec1, $vec2) = @_;

  my $dot = $$vec1[0] * $$vec2[0] + $$vec1[1] * $$vec2[1] + $$vec1[2] * $$vec2[2];
  return $dot;
  }

sub _vecLength
  {
  my $vec = shift @_;

  my $size = sqrt(($$vec[0])**2 + ($$vec[1])**2 + ($$vec[2])**2);

  return $size;
  }


sub _arccos
  {

  atan2( sqrt(1 - $_[0] * $_[0]), $_[0] )
  }

sub coordinates
  {
  my $self = shift @_;

  return join (", ", $self->{x}, $self->{y}, $self->{z})
  }

sub transform
  {
  my $self = shift @_;
  my ($matrix, $crystA, $crystB, $crystC) = @_;

  my $xnew = $$matrix[1][1] * $self->{x} + $$matrix[1][2] * $self->{y} + $$matrix[1][3] * $self->{z} + $$matrix[4][1] + $crystA;
  my $ynew = $$matrix[2][1] * $self->{x} + $$matrix[2][2] * $self->{y} + $$matrix[2][3] * $self->{z} + $$matrix[4][2] + $crystB;
  my $znew = $$matrix[3][1] * $self->{x} + $$matrix[3][2] * $self->{y} + $$matrix[3][3] * $self->{z} + $$matrix[4][3] + $crystC;

#print "$xnew, $ynew, $znew\n";
  my $newAtom = {%$self, "x" => $xnew, "y" => $ynew, "z" => $znew, "chainID" => join("", "#", $self->{chainID})} ; 
  return bless $newAtom, ref $self;
  }

sub transformBio
  {
  my $self = shift @_;
  my $matrix = shift @_;

  my $xnew = $$matrix[1][1] * $self->{x} + $$matrix[1][2] * $self->{y} + $$matrix[1][3] * $self->{z} + $$matrix[4][1];
  my $ynew = $$matrix[2][1] * $self->{x} + $$matrix[2][2] * $self->{y} + $$matrix[2][3] * $self->{z} + $$matrix[4][2];
  my $znew = $$matrix[3][1] * $self->{x} + $$matrix[3][2] * $self->{y} + $$matrix[3][3] * $self->{z} + $$matrix[4][3];

  my $newAtom = {%$self, "x" => $xnew, "y" => $ynew, "z" => $znew, "chainID" => join("", "*", $self->{chainID})} ;
  return bless $newAtom, ref $self;
  }

