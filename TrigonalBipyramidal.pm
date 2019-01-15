## TrigonalBipyramidal.pm 
############################################################################################
###
###   Written by Sen Yao, 07/20/2016
###   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
###
############################################################################################

package TrigonalBipyramidal;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 5
			  );

our $expectedAngle90 = 90;
our $expectedAngle120 = 120;
our $expectedAngle180 = 180;


our $invCorrM = [
	[ 0.802,-0.151,	-0.151,	 0.087,	-0.294,	-0.294,	 0.000,	 0.000,	 0.000],
	[-0.151, 0.802,	-0.151,	-0.294,	 0.087,	-0.294,	 0.000,	 0.000,	 0.000],
	[-0.151,-0.151,	 0.802,	-0.294,	-0.294,	 0.087,	 0.000,	 0.000,	 0.000],
	[ 0.087,-0.294,	-0.294,	 0.802,	-0.151,	-0.151,	 0.000,	 0.000,	 0.000],
	[-0.294, 0.087,	-0.294,	-0.151,	 0.802,	-0.151,	 0.000,	 0.000,	 0.000],
	[-0.294,-0.294,	 0.087,	-0.151,	-0.151,	 0.802,	 0.000,	 0.000,	 0.000],
	[ 0.000, 0.000,	 0.000,	 0.000,	 0.000,	 0.000,	 0.444,	-0.222,	-0.222],
	[ 0.000, 0.000,	 0.000,	 0.000,	 0.000,	 0.000,	-0.222,	 0.444,	-0.222],
	[ 0.000, 0.000,	 0.000,	 0.000,	 0.000,	 0.000,	-0.222,	-0.222,	 0.444]
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
    my @threeAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$twoAtoms)); } (@$combo);
    push (@$orderedCombos, [@$twoAtoms, @threeAtoms]);
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

  my ($mean90, $mean120, $mean180, $varianceOrN90, $varianceOrN120, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN120 = $self->numAngles() ;
    $varianceOrN180 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; #ideal mean, overal variance
    map { $angleTestStat += (($_ - $expectedAngle120) ** 2) / $varianceOrN120 ;} ($self->calcAllAngles120($combo)) ;
    map { $angleTestStat += (($_ - $expectedAngle180) ** 2) / $varianceOrN180 ;} ($self->calcAllAngles180($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean90 = ($$angleStats{"trigonalBipyramidal"}{"90"}{"mean"})? ($$angleStats{"trigonalBipyramidal"}{"90"}{"mean"}) : $expectedAngle90 ;
      $mean120 = ($$angleStats{"trigonalBipyramidal"}{"120"}{"mean"})? ($$angleStats{"trigonalBipyramidal"}{"120"}{"mean"}) : $expectedAngle120 ;
      #$varianceOrN90 = ($$angleStats{"trigonalBipyramidal"}{"90"}{"variance"})? ($$angleStats{"trigonalBipyramidal"}{"90"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      #$varianceOrN120 = ($$angleStats{"trigonalBipyramidal"}{"120"}{"variance"})? ($$angleStats{"trigonalBipyramidal"}{"120"}{"variance"}) : $$angleStats{"variance"}; 
      $varianceOrN90 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN120 = $$angleStats{"variance"};
      }
    else 
      {
      $varianceOrN90 =  $$angleStats{"trigonalBipyramidal"}{"90"}{"variance"};
      $varianceOrN120 = $$angleStats{"trigonalBipyramidal"}{"120"}{"variance"};
      }

    #my @means = ($expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle120, $expectedAngle120, $expectedAngle120); ## ideal mean
    my @means = ($mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean120, $mean120, $mean120); ## major mean
    my @angles = ($self->calcAllAngles90($combo), $self->calcAllAngles120($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std90 = 1/sqrt($varianceOrN90);
    my $std120 = 1/sqrt($varianceOrN120);
    my $invStds = [$std90, $std90, $std90, $std90, $std90, $std90, $std120, $std120, $std120];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    #print "mean 90, $expect90, mean 120, $expect120\nangles: ";
    #print map {"$_, "; } (@angles);
    #print "\nstd90, $std90, std 120, $std120\n";
    #print "covariance: ", $chiStat->element(1,1);#, "; previous: $angleTestStat\n";

    return $chiStat;
    }
  }


sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles90());
  push (@anglelist, $self->calcAllAngles120());
  push (@anglelist, $self->calcAllAngles180());

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"trigonalBipyramidal"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"trigonalBipyramidal"}{"120"}}, $_->calcAllAngles120());
  push (@{$$angleStats{"trigonalBipyramidal"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


sub calcAllAngles180
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4]];
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < $#$axial; $x++) # axial atom to axial atom angles, 180 as ideal
    {
    for(my $y = $x+1; $y < @$axial; $y++)
      {
      push (@angles, $center->angle($$axial[$x], $$axial[$y])) ;
      } 
    }

  return @angles;
  }

sub calcAllAngles120
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4]];
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < $#$planar; $x++) # planar atom to planar atom angles, 120 as ideal
    {
    for(my $y = $x+1; $y < @$planar; $y++)
      {
      push (@angles, $center->angle($$planar[$x], $$planar[$y]));
      }
    }

  return @angles;
  }

sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0], $$combo[1]];
  my $planar = [$$combo[2], $$combo[3], $$combo[4]];
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < @$axial; $x++) # axial atom to planar atom angles, 90 as ideal
    {
    for(my $y = 0; $y < @$planar; $y++)
      {
      push (@angles, $center->angle($$axial[$x], $$planar[$y]));
      }
    }

  return @angles;
  }


