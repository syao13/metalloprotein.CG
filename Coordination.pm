## Coordination.pm

package Coordination;
use strict;
use Atom;
use AtomShell;
use Distributions;
use Math::MatrixReal;

our @defaultDataMembers = (
                            "shellObj" => 0, # ref to AtomShell object
                            "combos" => 0    # array of possible atom combinations for coordination 
			    );

our $slope = 0.0562;

sub new
  {
  my $class = shift @_;
  my $self = { @defaultDataMembers, @_ };
  bless $self, ref $class || $class;

  $self->calculateCombos($self->{numAtoms});
  
  return $self;
  }

## calculateCombos
####  returns all nonredundant combinations of a desired size from residue list, with reasonable atom-atom distances
####    &restricted_combinations(coordination_object, size);
sub calculateCombos
  {
  my $self = shift @_;
  my $num = shift @_;
  my $upper_cutoff = (@_) ? shift @_ : 6.0;

  my $combos = [];
  if ($self->{numAtoms} <= @{$self->{shellObj}->{shell}})
    {
    my $n = 0; my $m = 0;
    foreach my $ligands (&_restrictedCombinations($num, @{$self->{shellObj}->{shell}}))
      {
      $m++;
      my @distances = &_distanceBetweenAtoms($ligands); #eliminates sets with unreasonable atom-atom distances
      if ( scalar(grep { $_ > 1.5 && $_ < $upper_cutoff; } (@distances)) == scalar @distances)
        { push @$combos, $ligands; $n++ }
      }
    }

  $self->{combos} = $combos; 
  }

## _restrictedCombinations
###  returns all nonredundant combinations of a desired size from a list
###    &restricted_combinations(size, list);
sub _restrictedCombinations
  {
  my $size = shift @_;
  return [] unless (@_ && $size);
  my $first = shift @_;
  my @rest1 = &_restrictedCombinations($size-1, @_);
  return (map { [$first, @$_] } (@rest1)) if ($size > @_);
  my @rest2 = &_restrictedCombinations($size, @_);

  return (@rest2, map { [$first, @$_] } (@rest1));
  }


sub distanceChi
  {
  my $self = shift @_;
  my $combo = shift @_;
  return 0 if (! @_);
  my $distanceStats = shift @_;

  my $distChi = 0;

  my $center = $self->{shellObj}->{center};
  for(my $x = 0; $x < @$combo; $x++)
    {
    my $element = $$combo[$x]->{element};
    my ($expect, $variance, $varianceOld);

    if (exists $$distanceStats{$element})
      {
      $expect = $$distanceStats{$element}{"mean"}; 
      $varianceOld = $$distanceStats{$element}{"variance"};

      my $resolution = ($$combo[$x]->{resolution} == -1)? 2.5 : $$combo[$x]->{resolution};
      my $adjStd = ($resolution - $$distanceStats{$element}{resolutionAvg}) * $slope + $$distanceStats{$element}{standardDeviation};
      $variance = $adjStd ** 2; 
      }
    else     
      {
      $expect = $$distanceStats{"average"}{"mean"};
      $varianceOld = $$distanceStats{"average"}{"variance"};
      my $resolution = ($$combo[$x]->{resolution} == -1)? 2.5 : $$combo[$x]->{resolution};
      my $adjStd = ($resolution - $$distanceStats{"average"}{resolutionAvg}) * $slope + $$distanceStats{"average"}{standardDeviation};
      $variance = $adjStd ** 2;
      }

    $distChi += ($center->distance($$combo[$x]) - $expect) ** 2 / $variance ;
    }

  return $distChi;
  }


sub bestDistChi
  {
  my $self = shift @_;
  my $angleDistStats = (@_) ? shift @_ : 0;

  my $bestStat;
  my $df = $self->{numAtoms} ;

  foreach my $combo ( @{$self->{combos}} )
    { 
    my $distChi = $self->distanceChi($combo, $$angleDistStats{"distance"}) ;
    my $probability = &Statistics::Distributions::chisqrprob($df, $distChi);

    if ( (! $bestStat) || $probability > $$bestStat{"probability"} )
      { $bestStat = { "ligands" => $combo, "probability" => $probability , "distChi" => $distChi }; }
    }

  $self->{bestCombo} = $bestStat;
  }

## If passing in $angleDistStats, it is calculating best Chi-square statistic, otherwise it is calculating best deviation 
sub bestTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $control = shift @_;
  my $threshold = shift @_;
  my $leaveOut = (@_) ? shift @_ : 0;
  my $angleDistStats = (@_) ? shift @_ : 0;

  my $bestStat;
  my $distChi = 0;
  my $df = 0;
  
  foreach my $combo ( @{$self->{combos}} )
    {
    #next if (grep {$_ < 68 ;} (&_anglesBetweenAtoms($combo, $self->{shellObj}->{center})) ); ## remove compressed angle combo. Used on 4.25 & 4.29
    next if ($control eq "c" && grep {$_ < $threshold;} (&_anglesBetweenAtoms($combo, $self->{shellObj}->{center})) ); ## remove compressed angle combo

    if ($angleDistStats) 
      {
      $distChi = $self->distanceChi($combo, $$angleDistStats{"distance"}) ;
      $df = $self->degreeFreedom() ;
      $df = $df - 1 if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut");
      }
    
    foreach my $orderedCombo (@{$self->orderedCombinations($combo)})
      {
      if ($type eq "dev")
	{
	my $testStatistics = $self->angleTestStatistic($type, $orderedCombo);

	if ( (! $bestStat) || $testStatistics < $$bestStat{"deviation"} )
	  { $bestStat = { "ligands" => $orderedCombo, "deviation" =>$testStatistics }; }
	}
      else
	{
	my $angleChi = $self->angleTestStatistic($type, $orderedCombo, $$angleDistStats{"angle"}, $leaveOut);
	
	my $testStatistics = $angleChi + $distChi;
	my $probability = &Statistics::Distributions::chisqrprob($df, $testStatistics);

	if ( (! $bestStat) || $probability > $$bestStat{"probability"} )
          { $bestStat = { "ligands" => $orderedCombo, "chiStatistic" => $testStatistics, "probability" => $probability , "angleChi" => $angleChi }; }
	}
      }
    }

  $self->{bestCombo} = $bestStat;
  }


sub numAngles
  {
  my $self = shift @_;

#  return ( $self->{numAtoms} * 2 - 3);
  return ($self->{numAtoms} * ($self->{numAtoms} - 1) / 2 ); 
  }


sub degreeFreedom
  { 
  my $self = shift @_;

  return ( $self->{numAtoms} * 3 - 3 );
  }


# distance_between_atoms
# #    calculates pairwise distances between a list of atoms
# # Parameters:
# #    $atoms - reference to array of atom info hash references
sub _distanceBetweenAtoms
  {
  my $atoms = shift @_;
  my @distances;
  for (my $x=0; $x < (scalar @$atoms) -1; $x++)
    {
    for (my $y=$x+1; $y < scalar @$atoms; $y++)
      {
      push @distances, $$atoms[$x]->distance($$atoms[$y]) ;
      }
    }
  return @distances;
  }

# angles_between_atoms
# # # #  calculates every atom to atom angle in a set of atoms
# # # #  &angles_between_atoms ( ref to array of atom information hashes, zinc information hash)
sub _anglesBetweenAtoms
  {
  my $atoms = shift @_;
  my $center = shift @_;
  my @angles;
  for (my $a=0; $a < (scalar @$atoms) -1; $a++)
    {
    for (my $b=$a+1; $b < scalar @$atoms; $b++)
      {
      push @angles, $center->angle($$atoms[$a],$$atoms[$b]);
      }
    }
  return @angles;
  }


sub covMatChi
  {
  my ($self, $diff, $invStds, $invCorr) = @_;

  my $diff1xN = Math::MatrixReal->new_from_rows([$diff]);
  my $diffNx1 = Math::MatrixReal->new_from_cols([$diff]);
  my $invStdM = Math::MatrixReal->new_diag($invStds);
  my $invCorrM = Math::MatrixReal->new_from_rows($invCorr);

  my $chi = ($diff1xN * ($invStdM * $invCorrM * $invStdM) * $diffNx1)->element(1,1);

  return $chi;
  }


sub bidentate
  {
  my $self = shift @_;
  return 0 if (! exists $self->{bestCombo});

  my %residue;
  foreach my $ligand ( @{$self->{bestCombo}->{ligands}})
    {
    my $key = $ligand->{residueName}.$ligand->resID();
    $residue{$key} += 1;
    }
  
  my $bi;
  foreach my $key (sort keys %residue)
    {
    if ($residue{$key} > 1)
      {
      if (length $bi == 0)
	{ $bi = $residue{$key}.".".(substr $key, 0, 3); }
      else
        { $bi = $bi.".".$residue{$key}.".".(substr $key, 0, 3) if (length $bi > 0); }
      }
    }

  return $bi if (length $bi > 0);
  #return "2.$_" if (grep {$residue{$_} == 2;} (keys %residue));
  return 1;
  }



sub ligandCombos
  {
  my $self = shift @_;
  my $comboLigands = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @ligands;
  foreach my $ligand (@$comboLigands)
    {
    push @ligands, $ligand->resID().".".$ligand->{residueName}.".".$ligand->{atomName};
    }

#  @ligands = sort @ligands;
  my (@chains, @residues, @atoms);
  foreach my $resAtom (@ligands)
    {
    my ($chain, $res, $atom) = split ('\.', $resAtom);
    push @chains, $chain;
    push @residues, $res;
    push @atoms, $atom;
    }
  
  my $amineN = (grep {$_ eq "N";} (@atoms))? 1: 0;
  return ( join (",", @chains), join (",", @residues), join (",", @atoms), $amineN);
  }


# for random forest
sub ligandAtomElement
  {
  my $self = shift @_;
  my $comboLigands = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @resAtomEle;
  my $center = $self->{shellObj}->{center};
 
  map {push @resAtomEle, $_->{residueName}.".".$_->{atomName}.".".$_->{element} ;} (@$comboLigands);
  map {push @resAtomEle, $center->distance($_) ;} (@$comboLigands);

#      push @resAtomEle, $$comboLigands[$x]->{residue_name}.".".$$comboLigands[$x]->{atom_name}.".".$$comboLigands[$x]->{element} ;
#      push @resAtomEle, $$comboLigands[$y]->{residue_name}.".".$$comboLigands[$y]->{atom_name}.".".$$comboLigands[$x]->{element} ;

  for (my $x = 0; $x < $#$comboLigands; $x++)
    {
    for (my $y = $x+1; $y < @$comboLigands; $y++)
      {
      if ($$comboLigands[$x]->resID() eq $$comboLigands[$y]->resID())
	{ push @resAtomEle, 1;}
      else
        { push @resAtomEle, 0;}
      }
    }

  return @resAtomEle;

  }


## Order angle as largest, middles, and opposite, with the ligands ordered as simply the two composing the largest and then the two composing the opposite
sub orderedAngles
  {
  my $self = shift @_;
  return 0 if (! exists $self->{bestCombo});
  my $combo = $self->{bestCombo}->{ligands};

  my @anglelist = $self->allAngles(); ## Angles ordered by the ligand position in combo
  @anglelist = sort {$b <=> $a} (@anglelist);

  my $center = $self->{shellObj}->{center};
  my %largest;
  FORLOOP: for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      my $angle = $center->angle($$combo[$x], $$combo[$y]);

      if ($angle == $anglelist[0])
        {
        $largest{$x} = 1;
        $largest{$y} = 1;
        last FORLOOP;
        }
      }
    }

  my @otherinds = grep { $largest{$_} != 1;} (0..(@$combo-1)); ## find the smallest opposite
  my $otherLigs = [map {$$combo[$_]} (grep {my $ind = $_; $largest{$ind} != 1;} (0..(@$combo-1)))];
  my @otherAngles = $self->allAngles($otherLigs);
  @otherAngles = sort {$a <=> $b} (@otherAngles);

  my %opposite;
  for (my $x = 2; $x < @otherAngles - 1 ; $x++)
    {
    for(my $y = $x+1; $y < @otherAngles; $y++)
      {
      my $angle = $center->angle($$combo[$otherinds[$x]], $$combo[$otherinds[$y]]);

      if ($angle == $otherAngles[0])
	{
        $opposite{$otherinds[$x]} = 1;
        $opposite{$otherinds[$y]} = 1;
	}
      }
    }

#print "@anglelist\n";
#print %largest, "\n";
#print "otherinds, @otherinds\n";
#print "otherligs, ", @$otherLigs, "\n";
#print "@otherAngles\n";
#print %opposite, "\n";

  my $newCombo = [];
  push @$newCombo, (map {$$combo[$_]} (keys %largest));
  push @$newCombo, (map {$$combo[$_]} (grep {my $ind = $_; $largest{$ind} != 1 && $opposite{$ind} != 1;} (0..(@$combo-1))));
  push @$newCombo, (map {$$combo[$_]} (keys %opposite));
  $self->{bestCombo}->{ligands} = $newCombo;

  return $self->allAngles();
  }


sub allAngles
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



















