## PentagonalBipyramidal.pm 

package PentagonalBipyramidal;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 7
			  );

our $expectedAngle72 = 72;
our $expectedAngle90 = 90;
our $expectedAngle144 = 144;
our $expectedAngle180 = 180;


our $invCorrM = [
		];

sub new
  {
  my $class = shift @_;
  my $self = Coordination->new(@defaultDataMembers, @_);

  bless $self, ref $class || $class;

  return $self;
  }


sub orderedCombinations
  {
  my $self = shift @_;
  my $combo = shift @_;

  my $orderedCombos;
  foreach my $twoAtoms (&Coordination::_restrictedCombinations(2, @$combo))
    {
    my @fiveAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$twoAtoms)); } (@$combo);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[1],$fiveAtoms[2],$fiveAtoms[3], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[1],$fiveAtoms[3],$fiveAtoms[2], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[1],$fiveAtoms[4],$fiveAtoms[2], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[2],$fiveAtoms[3],$fiveAtoms[1], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[2],$fiveAtoms[4],$fiveAtoms[1], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[3],$fiveAtoms[4],$fiveAtoms[1], $fiveAtoms[2]]);
    }

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean90, $mean180, $varianceOrN72, $varianceOrN90, $varianceOrN144, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN72 = $self->numAngles() ;
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN144 = $self->numAngles() ;
    $varianceOrN180 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle72) ** 2) / $varianceOrN72 ;} ($self->calcAllAngles72($combo)) ; #ideal mean, number of variable
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle144) ** 2) / $varianceOrN144 ;} ($self->calcAllAngles144($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle180) ** 2) / $varianceOrN180 ;} ($self->calcAllAngles180($combo)) ;

    return $angleTestStat;
    }
  else
    {
    print STDERR "ERROR, should not be here in bootstrap!\n";
    }
  }


sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < @$axial; $x++) # axial atom to planar atom angles, 90 is ideal
    {
    for(my $y = 0; $y < @$planar; $y++)
      {
      push (@angles, $center->angle($$axial[$x], $$planar[$y]));
      }
    }

  return @angles;
  }

sub calcAllAngles180
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$axial[0], $$axial[1])) ; # axial atom to axial atom angles, 180 as ideal

  return @angles;
  }

sub calcAllAngles72
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[1])); ## planar atom to adjacent planar atom angles, 72 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[2])); 
  push (@angles, $center->angle($$planar[2], $$planar[3])); 
  push (@angles, $center->angle($$planar[3], $$planar[4])); 
  push (@angles, $center->angle($$planar[4], $$planar[0])); 

  return @angles;
  }

sub calcAllAngles144
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[2])); ## planar atom to every other planar atom angles, 144 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[3])); 
  push (@angles, $center->angle($$planar[2], $$planar[4])); 
  push (@angles, $center->angle($$planar[3], $$planar[0])); 
  push (@angles, $center->angle($$planar[4], $$planar[1]));

  return @angles;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"pentagonalBipyramidal"}{"72"}}, $_->calcAllAngles72());
  push (@{$$angleStats{"pentagonalBipyramidal"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"pentagonalBipyramidal"}{"144"}}, $_->calcAllAngles144());
  push (@{$$angleStats{"pentagonalBipyramidal"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


