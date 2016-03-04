## TrigonalPrismaticV.pm 

package TrigonalPrismaticV;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 5
			  );

our $expectedAngle90 = 90;
our $expectedAngle70 = 70.6;
our $expectedAngle131 = 131.8;


our $invCorrM = [ 
	[ 0.655,  0.076,  0.076, -0.260, -0.082, -0.082,  0.164,  0.164,  0.045,  0.045],
	[ 0.076,  0.640,  0.307, -0.198,  0.174,  0.066,  0.185, -0.046, -0.250, -0.101],
	[ 0.076,  0.307,  0.640, -0.198,  0.066,  0.174, -0.046,  0.185, -0.101, -0.250],
	[-0.260, -0.198, -0.198,  0.486, -0.167, -0.167, -0.061, -0.061, -0.071, -0.071],
	[-0.082,  0.174,  0.066, -0.167,  0.665, -0.181,  0.013,  0.264, -0.010, -0.148],
	[-0.082,  0.066,  0.174, -0.167, -0.181,  0.665,  0.264,  0.013, -0.148, -0.010],
	[ 0.164,  0.185, -0.046, -0.061,  0.013,  0.264,  0.550, -0.058, -0.049, -0.186],
	[ 0.164, -0.046,  0.185, -0.061,  0.264,  0.013, -0.058,  0.550, -0.186, -0.049],
	[ 0.045, -0.250, -0.101, -0.071, -0.010, -0.148, -0.049, -0.186,  0.380, -0.044],
	[ 0.045, -0.101, -0.250, -0.071, -0.148, -0.010, -0.186, -0.049, -0.044,  0.380]
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
  foreach my $topAtoms (&Coordination::_restrictedCombinations(2, @$combo))
    {
    my @bottomAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$topAtoms)); } (@$combo);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[0], $bottomAtoms[1], $bottomAtoms[2]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[1], $bottomAtoms[2], $bottomAtoms[0]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[2], $bottomAtoms[0], $bottomAtoms[1]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[0], $bottomAtoms[2], $bottomAtoms[1]]);
    push (@$orderedCombos, [@$topAtoms, $bottomAtoms[1], $bottomAtoms[0], $bottomAtoms[2]]);
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
      $mean90 = ($$angleStats{"trigonalPrismaticVancancy"}{"90"}{"mean"})? ($$angleStats{"trigonalPrismaticVancancy"}{"90"}{"mean"}) : $expectedAngle90 ;
      $mean70 = ($$angleStats{"trigonalPrismaticVancancy"}{"70"}{"mean"})? ($$angleStats{"trigonalPrismaticVancancy"}{"70"}{"mean"}) : $expectedAngle70 ;
      $mean131 = ($$angleStats{"trigonalPrismaticVancancy"}{"131"}{"mean"})? ($$angleStats{"trigonalPrismaticVancancy"}{"131"}{"mean"}) : $expectedAngle131 ;
      $varianceOrN90 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN70 = $$angleStats{"variance"};
      $varianceOrN131 = $$angleStats{"variance"};
      }
    else 
      {
      $varianceOrN90 =  $$angleStats{"trigonalPrismaticVancancy"}{"90"}{"variance"};
      $varianceOrN70 = $$angleStats{"trigonalPrismaticVancancy"}{"70"}{"variance"};
      }

    my @means = ($mean90, $mean90, $mean90, $mean90, $mean70, $mean70, $mean131, $mean131, $mean131, $mean131); ## major mean
    my @angles = ($self->calcAllAngles90($combo), $self->calcAllAngles70($combo), $self->calcAllAngles131($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    my $std90 = 1/sqrt($varianceOrN90);
    my $std70 = 1/sqrt($varianceOrN70);
    my $std131 = 1/sqrt($varianceOrN131);
    my $invStds = [$std90, $std90, $std90, $std90, $std70, $std70, $std131, $std131, $std131, $std131];

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

  push (@{$$angleStats{"trigonalPrismaticVancancy"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"trigonalPrismaticVancancy"}{"70"}}, $_->calcAllAngles70());
  push (@{$$angleStats{"trigonalPrismaticVancancy"}{"131"}}, $_->calcAllAngles131());

  return 0;
  }


sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[2], $$combo[3])) ;
  push (@angles, $center->angle($$combo[3], $$combo[4])) ;
  push (@angles, $center->angle($$combo[4], $$combo[5])) ;
  push (@angles, $center->angle($$combo[0], $$combo[1])) ;

  return @angles;
  }

sub calcAllAngles70
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[2])) ;
  push (@angles, $center->angle($$combo[1], $$combo[3])) ;

  return @angles;
  }

sub calcAllAngles131
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[1], $$combo[2])) ;
  push (@angles, $center->angle($$combo[0], $$combo[3])) ;
  push (@angles, $center->angle($$combo[1], $$combo[4])) ;
  push (@angles, $center->angle($$combo[0], $$combo[4])) ;

  return @angles;
  }


