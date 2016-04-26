## TrigonalPrismatic.pm 

package TrigonalPrismatic;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 6
			  );

our $expectedAngle90 = 90;
our $expectedAngle70 = 70.6;
our $expectedAngle131 = 131.8;


our $invCorrM = [ 
	[ 0.419,  0.093,  0.093, -0.176, -0.110, -0.110, -0.104, -0.104,  0.061,  0.031, -0.127,  0.031, -0.127,  0.049,  0.049],
	[ 0.093,  0.419,  0.093, -0.110, -0.176, -0.110,  0.061, -0.104, -0.104,  0.049,  0.049, -0.127,  0.031, -0.127,  0.031],
	[ 0.093,  0.093,  0.419, -0.110, -0.110, -0.176, -0.104,  0.061, -0.104, -0.127,  0.031,  0.049,  0.049,  0.031, -0.127],
	[-0.176, -0.110, -0.110,  0.419,  0.093,  0.093, -0.104, -0.104,  0.061,  0.031,  0.049,  0.031,  0.049, -0.127, -0.127],
	[-0.110, -0.176, -0.110,  0.093,  0.419,  0.093,  0.061, -0.104, -0.104, -0.127, -0.127,  0.049,  0.031,  0.049,  0.031],
	[-0.110, -0.110, -0.176,  0.093,  0.093,  0.419, -0.104,  0.061, -0.104,  0.049,  0.031, -0.127, -0.127,  0.031,  0.049],
	[-0.104,  0.061, -0.104, -0.104,  0.061, -0.104,  0.478, -0.135, -0.135,  0.077,  0.077,  0.077, -0.088,  0.077, -0.088],
	[-0.104, -0.104,  0.061, -0.104, -0.104,  0.061, -0.135,  0.478, -0.135,  0.077, -0.088,  0.077,  0.077, -0.088,  0.077],
	[ 0.061, -0.104, -0.104,  0.061, -0.104, -0.104, -0.135, -0.135,  0.478, -0.088,  0.077, -0.088,  0.077,  0.077,  0.077],
	[ 0.031,  0.049, -0.127,  0.031, -0.127,  0.049,  0.077,  0.077, -0.088,  0.311, -0.113, -0.034, -0.005, -0.005, -0.113],
	[-0.127,  0.049,  0.031,  0.049, -0.127,  0.031,  0.077, -0.088,  0.077, -0.113,  0.311, -0.005, -0.113, -0.034, -0.005],
	[ 0.031, -0.127,  0.049,  0.031,  0.049, -0.127,  0.077,  0.077, -0.088, -0.034, -0.005,  0.311, -0.113, -0.113, -0.005],
	[-0.127,  0.031,  0.049,  0.049,  0.031, -0.127, -0.088,  0.077,  0.077, -0.005, -0.113, -0.113,  0.311, -0.005, -0.034],
	[ 0.049, -0.127,  0.031, -0.127,  0.049,  0.031,  0.077, -0.088,  0.077, -0.005, -0.034, -0.113, -0.005,  0.311, -0.113],
	[ 0.049,  0.031, -0.127, -0.127,  0.031,  0.049, -0.088,  0.077,  0.077, -0.113, -0.005, -0.005, -0.034, -0.113,  0.311]
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
  my %tops;
  foreach my $topAtoms (&Coordination::_restrictedCombinations(3, @$combo))
    {
    next if ($tops{join ("", sort (@$topAtoms))});

    my @bottomAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$topAtoms)); } (@$combo);
    my $bottomString = join ("", sort (@bottomAtoms));
    $tops{$bottomString} = 1;

    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[0], $bottomAtoms[1], $bottomAtoms[2]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[0], $bottomAtoms[2], $bottomAtoms[1]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[1], $bottomAtoms[0], $bottomAtoms[2]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[1], $bottomAtoms[2], $bottomAtoms[0]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[2], $bottomAtoms[0], $bottomAtoms[1]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[2], $bottomAtoms[1], $bottomAtoms[0]]);
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

  my ($mean90, $mean70, $mean131, $varianceOrN90, $varianceOrN70, $varianceOrN131);

  if ($type eq "dev")
    {
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN70 = $self->numAngles() ;
    $varianceOrN131 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; #ideal mean, overal variance
    map { $angleTestStat += (($_ - $expectedAngle70) ** 2) / $varianceOrN70 ;} ($self->calcAllAngles70($combo)) ;
    map { $angleTestStat += (($_ - $expectedAngle131) ** 2) / $varianceOrN131 ;} ($self->calcAllAngles131($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean90 = ($$angleStats{"trigonalPrismatic"}{"90"}{"mean"})? ($$angleStats{"trigonalPrismatic"}{"90"}{"mean"}) : $expectedAngle90 ;
      $mean70 = ($$angleStats{"trigonalPrismatic"}{"70"}{"mean"})? ($$angleStats{"trigonalPrismatic"}{"70"}{"mean"}) : $expectedAngle70 ;
      $mean131 = ($$angleStats{"trigonalPrismatic"}{"131"}{"mean"})? ($$angleStats{"trigonalPrismatic"}{"131"}{"mean"}) : $expectedAngle131 ;
      $varianceOrN90 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN70 = $$angleStats{"variance"};
      $varianceOrN131 = $$angleStats{"variance"};
      }
    else 
      {
      $varianceOrN90 =  $$angleStats{"trigonalPrismatic"}{"90"}{"variance"};
      $varianceOrN70 = $$angleStats{"trigonalPrismatic"}{"70"}{"variance"};
      }

    #my @means = ($expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle70, $expectedAngle70, $expectedAngle70); ## ideal mean
    my @means = ($mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean70, $mean70, $mean70, $mean131, $mean131, $mean131, $mean131, $mean131, $mean131); ## major mean
    my @angles = ($self->calcAllAngles90($combo), $self->calcAllAngles70($combo), $self->calcAllAngles131($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std90 = 1/sqrt($varianceOrN90);
    my $std70 = 1/sqrt($varianceOrN70);
    my $std131 = 1/sqrt($varianceOrN131);
    my $invStds = [$std90, $std90, $std90, $std90, $std90, $std90, $std70, $std70, $std70, $std131, $std131, $std131, $std131, $std131, $std131];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    #print "mean 90, $expect90, mean 70, $expect70\nangles: ";
    #print map {"$_, "; } (@angles);
    #print "\nstd90, $std90, std 70, $std70\n";
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
  push (@anglelist, $self->calcAllAngles70());
  push (@anglelist, $self->calcAllAngles131());

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"trigonalPrismatic"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"trigonalPrismatic"}{"70"}}, $_->calcAllAngles70());
  push (@{$$angleStats{"trigonalPrismatic"}{"131"}}, $_->calcAllAngles131());

  return 0;
  }


sub calcAllAngles131
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[4])) ;
  push (@angles, $center->angle($$combo[0], $$combo[5])) ;
  push (@angles, $center->angle($$combo[1], $$combo[3])) ;
  push (@angles, $center->angle($$combo[1], $$combo[5])) ;
  push (@angles, $center->angle($$combo[2], $$combo[3])) ;
  push (@angles, $center->angle($$combo[2], $$combo[4])) ;

  return @angles;
  }

sub calcAllAngles70
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[3])) ;
  push (@angles, $center->angle($$combo[1], $$combo[4])) ;
  push (@angles, $center->angle($$combo[2], $$combo[5])) ;

  return @angles;
  }

sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[1])) ;
  push (@angles, $center->angle($$combo[1], $$combo[2])) ;
  push (@angles, $center->angle($$combo[2], $$combo[0])) ;
  push (@angles, $center->angle($$combo[3], $$combo[4])) ;
  push (@angles, $center->angle($$combo[4], $$combo[5])) ;
  push (@angles, $center->angle($$combo[5], $$combo[3])) ;

  return @angles;
  }


