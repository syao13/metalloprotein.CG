## HexagonalBipyramidalVP.pm

package HexagonalBipyramidalVP;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 7
			  );

our $expectedAngle60 = 60;
our $expectedAngle90 = 90;
our $expectedAngle120 = 120;
our $expectedAngle180 = 180;


our $invCorrM = [
	[ 0.302, -0.151,  0.071, -0.142,  0.151, -0.080, -0.071, -0.080,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.151,  0.231, -0.080,  0.071,  0.080,  0.151, -0.009, -0.071,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.071, -0.080,  0.231, -0.151, -0.009,  0.151,  0.080, -0.071,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.142,  0.071, -0.151,  0.302, -0.071, -0.080,  0.151, -0.080,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.151,  0.080, -0.009, -0.071,  0.231,  0.071, -0.080, -0.151,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.080,  0.151,  0.151, -0.080,  0.071,  0.302,  0.071, -0.142,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.071, -0.009,  0.080,  0.151, -0.080,  0.071,  0.231, -0.151,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.080, -0.071, -0.071, -0.080, -0.151, -0.142, -0.151,  0.302,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.540, -0.023, -0.062, -0.040,  0.023, -0.248,  0.065, -0.188, -0.252, -0.065],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.023,  0.540,  0.063,  0.023, -0.040,  0.065, -0.248,  0.187, -0.065, -0.252],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.062,  0.063,  0.625,  0.063, -0.063, -0.188,  0.187, -0.125,  0.188, -0.188],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.040,  0.023,  0.063,  0.540, -0.023, -0.252, -0.065,  0.187, -0.248,  0.065],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.023, -0.040, -0.062, -0.023,  0.540, -0.065, -0.252, -0.188,  0.065, -0.248],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.248,  0.065, -0.188, -0.252, -0.065,  0.540, -0.023, -0.062, -0.040,  0.023],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.065, -0.248,  0.187, -0.065, -0.252, -0.023,  0.540,  0.063,  0.023, -0.040],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.187,  0.187, -0.125,  0.188, -0.188, -0.063,  0.063,  0.625,  0.062, -0.063],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.252, -0.065,  0.187, -0.248,  0.065, -0.040,  0.023,  0.062,  0.540, -0.023],
	[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, -0.065, -0.252, -0.187,  0.065, -0.248,  0.023, -0.040, -0.062, -0.023,  0.540]
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

    foreach my $idx1 (1..4)
      {
      foreach my $idx2 (1..4)
        {
        next if ($idx1 == $idx2);
        foreach my $idx3 (1..4)
          {
          next if ($idx1 == $idx3 || $idx2 == $idx3);
          my @rest = grep {$_ != $idx1 && $_ != $idx2 && $_ != $idx3} (1..4);

          push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[0], $fiveAtoms[$idx1], $fiveAtoms[$idx2], $fiveAtoms[$idx3], $fiveAtoms[$rest[0]]]);
          push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[$idx1], $fiveAtoms[0], $fiveAtoms[$idx2], $fiveAtoms[$idx3], $fiveAtoms[$rest[0]]]);
          }
        }
      }
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[2], $fiveAtoms[0], $fiveAtoms[3], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[2], $fiveAtoms[0], $fiveAtoms[4], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[3], $fiveAtoms[0], $fiveAtoms[2], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[3], $fiveAtoms[0], $fiveAtoms[4], $fiveAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[4], $fiveAtoms[0], $fiveAtoms[2], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[1], $fiveAtoms[4], $fiveAtoms[0], $fiveAtoms[3], $fiveAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[2], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[3], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[2], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[4], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[3], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[2], $fiveAtoms[4]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[3], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[4], $fiveAtoms[2]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[4], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[2], $fiveAtoms[3]]);
    push (@$orderedCombos, [@$twoAtoms, $fiveAtoms[4], $fiveAtoms[1], $fiveAtoms[0], $fiveAtoms[3], $fiveAtoms[2]]);
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

  my ($mean60, $mean120, $mean90, $mean180, $varianceOrN60, $varianceOrN90, $varianceOrN120, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN60 = $self->numAngles() ;
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN120 = $self->numAngles() ;
    $varianceOrN180 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle60) ** 2) / $varianceOrN60 ;} ($self->calcAllAngles60($combo)) ; #ideal mean, number of variable
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle120) ** 2) / $varianceOrN120 ;} ($self->calcAllAngles120($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle180) ** 2) / $varianceOrN180 ;} ($self->calcAllAngles180($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean60 = ($$angleStats{"HexagonalBipyramidal"}{"60"}{"mean"})? ($$angleStats{"HexagonalBipyramidal"}{"60"}{"mean"}) : $expectedAngle60; ## major coordinaiton mean
      $mean90 = ($$angleStats{"HexagonalBipyramidal"}{"90"}{"mean"})? ($$angleStats{"HexagonalBipyramidal"}{"90"}{"mean"}) : $expectedAngle90;
      $mean120 = ($$angleStats{"HexagonalBipyramidal"}{"120"}{"mean"})? ($$angleStats{"HexagonalBipyramidal"}{"120"}{"mean"}) : $expectedAngle120;
      $varianceOrN60 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN90 = $$angleStats{"variance"}; 
      $varianceOrN120 = $$angleStats{"variance"}; 
      }
    else
      {
      $varianceOrN60 =  $$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"60"}{"variance"};
      $varianceOrN90 =  $$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"90"}{"variance"};
      $varianceOrN120 =  $$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"120"}{"variance"};
      }

    my @means = ($mean60, $mean60, $mean60, $mean60, $mean120, $mean120, $mean120, $mean120, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90); 
    my @angles = ($self->calcAllAngles60($combo), $self->calcAllAngles120($combo), $self->calcAllAngles90($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std60 = 1/sqrt($varianceOrN60);
    my $std90 = 1/sqrt($varianceOrN90);
    my $std120 = 1/sqrt($varianceOrN120);
    my $invStds = [$std60, $std60, $std60, $std60, $std120, $std120, $std120, $std120, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90];

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
  push (@angles, $center->angle($$planar[0], $$planar[3])) ; # axial atom to axial atom angles, 180 as ideal
  push (@angles, $center->angle($$planar[1], $$planar[4])) ; # axial atom to axial atom angles, 180 as ideal

  return @angles;
  }

sub calcAllAngles60
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[1])); ## planar atom to adjacent planar atom angles, 60 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[2])); 
  push (@angles, $center->angle($$planar[2], $$planar[3])); 
  push (@angles, $center->angle($$planar[3], $$planar[4])); 

  return @angles;
  }

sub calcAllAngles120
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4], $$combo[5], $$combo[6]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[2])); ## planar atom to every other planar atom angles, 120 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[3])); 
  push (@angles, $center->angle($$planar[2], $$planar[4])); 
  push (@angles, $center->angle($$planar[4], $$planar[0]));

  return @angles;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"60"}}, $_->calcAllAngles60());
  push (@{$$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"120"}}, $_->calcAllAngles120());
  push (@{$$angleStats{"HexagonalBipyramidalVacancyPlanar"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


