## TrigonalBipyramidalVA.pm 

package TrigonalBipyramidalVA;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 4 
			  );

our $expectedAngle90 = 90;
our $expectedAngle120 = 120;

our $invCorrM = [
	[1.2,	0.4,	0.4,	0.0000000,	0.0000000,	0.0000000],
	[0.4,	1.2,	0.4,	0.0000000,	0.0000000,	0.0000000],
	[0.4,	0.4,	1.2,	0.0000000,	0.0000000,	0.0000000],
	[0.0,	0.0,	0.0,	0.4444444,	-0.2222222,	-0.2222222],
	[0.0,	0.0,	0.0,	-0.2222222,	0.4444444,	-0.2222222],
	[0.0,	0.0,	0.0,	-0.2222222,	-0.2222222,	0.4444444]
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
  push (@$orderedCombos, [$$combo[0], $$combo[1], $$combo[2], $$combo[3]]);
  push (@$orderedCombos, [$$combo[1], $$combo[0], $$combo[2], $$combo[3]]);
  push (@$orderedCombos, [$$combo[2], $$combo[0], $$combo[1], $$combo[3]]);
  push (@$orderedCombos, [$$combo[3], $$combo[0], $$combo[1], $$combo[2]]);

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean90, $mean120, $varianceOrN90, $varianceOrN120);

  if ($type eq "dev")
    {
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN120 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; #ideal mean, overal variance
    map { $angleTestStat += (($_ - $expectedAngle120) ** 2) / $varianceOrN120 ;} ($self->calcAllAngles120($combo)) ;

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
    elsif ($type eq "ownStats")
      {
      $mean90 = $$angleStats{"trigonalBipyramidalVacancyAxial"}{"90"}{"mean"} ;
      $mean120 = $$angleStats{"trigonalBipyramidalVacancyAxial"}{"120"}{"mean"} ;
      $varianceOrN90 =  $$angleStats{"trigonalBipyramidalVacancyAxial"}{"90"}{"variance"};
      $varianceOrN120 = $$angleStats{"trigonalBipyramidalVacancyAxial"}{"120"}{"variance"};
      }

    #my @means = ($expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle120, $expectedAngle120, $expectedAngle120); ## ideal mean
    my @means = ($mean90, $mean90, $mean90, $mean120, $mean120, $mean120); ## major mean
    my @angles = ($self->calcAllAngles90($combo), $self->calcAllAngles120($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    ## set the smallest angle's diff to zero
    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        {if ($angles[$x] < $angles[$n]) {$n = $x;} } 
      $$diff[$n] = 0;
      }

    my $std90 = 1/sqrt($varianceOrN90);
    my $std120 = 1/sqrt($varianceOrN120);
    my $invStds = [$std90, $std90, $std90, $std120, $std120, $std120];
    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    return $chiStat;
    }
  }



sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist ;
  push (@anglelist, $self->calcAllAngles90());
  push (@anglelist, $self->calcAllAngles120());

  my $pair1 = [sort {$a <=> $b} ($anglelist[0], $anglelist[5])];
  my $pair2 = [sort {$a <=> $b} ($anglelist[1], $anglelist[4])];
  my $pair3 = [sort {$a <=> $b} ($anglelist[2], $anglelist[3])];

  my @pairlist = sort { $$a[0] <=> $$b[0]; } ($pair1, $pair2, $pair3);
  @anglelist = ( @{$pairlist[0]}, @{$pairlist[1]}, @{$pairlist[2]} );

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"trigonalBipyramidalVacancyAxial"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"trigonalBipyramidalVacancyAxial"}{"120"}}, $_->calcAllAngles120());

  return 0;
  }


sub calcAllAngles120
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $planar = [$$combo[1], $$combo[2], $$combo[3]];
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
  my $axial = $$combo[0];
  my $planar = [$$combo[1], $$combo[2], $$combo[3]];
  my $center = $self->{shellObj}->{center};

  # axial atom to planar atom angles, 90 as ideal
  for(my $x = 0; $x < @$planar; $x++)
    {
    push (@angles, $center->angle($axial, $$planar[$x]));
    }
  
  return @angles;
  }


sub calcAllAngles
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my $center = $self->{shellObj}->{center};

  my @angles;
  for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      push (@angles, $center->angle($$combo[$x], $$combo[$y]));
      }
    }

  return @angles;
  }


