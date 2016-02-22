## PentagonalBipyramidal.pm 

package PentagonalBipyramidalVP;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 6
			  );

our $expectedAngle72 = 72;
our $expectedAngle90 = 90;
our $expectedAngle144 = 144;
our $expectedAngle180 = 180;


our $invCorrM = [
	[ 0.250, -0.125,  0.000,  0.125, -0.125, -0.125,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.125,  0.250, -0.125,  0.125,  0.125,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.000, -0.125,  0.250, -0.125,  0.125, -0.125,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.125,  0.125, -0.125,  0.250,  0.000, -0.125,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.125,  0.125,  0.125,  0.000,  0.250, -0.125,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.125,  0.000, -0.125, -0.125, -0.125,  0.250,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.574, -0.048, -0.106, -0.016, -0.195,  0.003, -0.299, -0.189],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.048,  0.714,  0.183, -0.106,  0.003, -0.020,  0.289, -0.299],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.106,  0.183,  0.714, -0.048, -0.299,  0.289, -0.020,  0.003],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.016, -0.106, -0.048,  0.574, -0.189, -0.299,  0.003, -0.195],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.195,  0.003, -0.299, -0.189,  0.574, -0.048, -0.106, -0.016],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.003, -0.020,  0.289, -0.299, -0.048,  0.714,  0.183, -0.106],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.299,  0.289, -0.020,  0.003, -0.106,  0.183,  0.714, -0.048],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.189, -0.299,  0.003, -0.195, -0.016, -0.106, -0.048,  0.574]
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
    my @fourAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$twoAtoms)); } (@$combo);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[1],$fourAtoms[2],$fourAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[1],$fourAtoms[3],$fourAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[2],$fourAtoms[1],$fourAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[2],$fourAtoms[3],$fourAtoms[1]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[3],$fourAtoms[1],$fourAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[0], $fourAtoms[3],$fourAtoms[2],$fourAtoms[1]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[1], $fourAtoms[0],$fourAtoms[3],$fourAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[1], $fourAtoms[0],$fourAtoms[2],$fourAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[1], $fourAtoms[2],$fourAtoms[0],$fourAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[1], $fourAtoms[3],$fourAtoms[0],$fourAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[2], $fourAtoms[0],$fourAtoms[1],$fourAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fourAtoms[2], $fourAtoms[1],$fourAtoms[0],$fourAtoms[3]]);
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

  my ($mean72, $mean144, $mean90, $mean180, $varianceOrN72, $varianceOrN90, $varianceOrN144, $varianceOrN180);

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
    if ($type eq "chi")
      {
      $mean72 = ($$angleStats{"pentagonalBipyramidal"}{"72"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"72"}{"mean"}) : $expectedAngle90; ## major coordinaiton mean
      $mean90 = ($$angleStats{"pentagonalBipyramidal"}{"90"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"90"}{"mean"}) : $expectedAngle90;
      $mean144 = ($$angleStats{"pentagonalBipyramidal"}{"144"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"144"}{"mean"}) : $expectedAngle90;
      $varianceOrN72 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN90 = $$angleStats{"variance"}; 
      $varianceOrN144 = $$angleStats{"variance"}; 
      }
    else
      {
      $varianceOrN72 =  $$angleStats{"pentagonalBipyramidal"}{"72"}{"variance"};
      $varianceOrN90 =  $$angleStats{"pentagonalBipyramidal"}{"90"}{"variance"};
      $varianceOrN144 =  $$angleStats{"pentagonalBipyramidal"}{"144"}{"variance"};
      }

    my @means = ($mean72, $mean72, $mean72, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean144, $mean144, $mean144); 
    my @angles = ($self->calcAllAngles72($combo), $self->calcAllAngles90($combo), $self->calcAllAngles144($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    my $std72 = 1/sqrt($varianceOrN72);
    my $std90 = 1/sqrt($varianceOrN90);
    my $std144 = 1/sqrt($varianceOrN144);
    my $invStds = [$std72, $std72, $std72, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std144, $std144, $std144];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    #print "mean 90, $expect90\nangles: ";
    ##print map {"$_, "; } (@angles);
    ##print "\nstd90, $std90\n";
    ##print "covariance: ", $chiStat;#, "; previous: $angleTestStat\n";
    
    return $chiStat;
    }
  }


sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5]];
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
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5]];
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
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[1])); ## planar atom to adjacent planar atom angles, 72 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[2])); 
  push (@angles, $center->angle($$planar[2], $$planar[3])); 

  return @angles;
  }

sub calcAllAngles144
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[2])); ## planar atom to every other planar atom angles, 144 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[3])); 
  push (@angles, $center->angle($$planar[3], $$planar[0])); 

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


