##!/usr/bin/perl 
#===============================================================================
#
#         FILE:  ZnCGanalysis.pm
#
#        USAGE:  ./ZnCGanalysis.pm 
#
#  DESCRIPTION:  A comprehensive process object that makes uses of lots of the
#  		 other objects and conducts the (zinc) coordination classification 
#  		 task
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Sen Yao 
#      COMPANY:  
#      VERSION:  1.0 
#      CREATED:  4/2/2013 02:07:30 PM
#     REVISION:  ---
#===============================================================================

package ZnCGanalysis;
use strict;
use BasicTable2 qw(:ALL);
use PDBEntry;
use Atom;
use AtomShell;
use Residue;
use Coordination;
use RawStatistics;
use Sequence;

use Tetrahedral;
use TrigonalBipyramidal;
use Octahedral;
use TrigonalPlanar;
use TetrahedralV;
use TrigonalBipyramidalVA;
use TrigonalBipyramidalVP;
use SquarePlanar;
use SquarePyramidalV;
use SquarePyramidal;

use Clone 'clone';

#use Data::Dumper::Concise; # human readable, code in iaCoordination
#use JSON; # For other programs to read back, code in iaCoordination 
#use JSON -convertBlessedUniversally;

our @defaultDataMembers = (
                          "pathsFile" => 0,  # a file that contains pdb-file paths
			  "element" => 0,
			  "shells" => 0,  # ref to AtomShell object
			  "coordinations" => 0,  
			  "stats" => 0,
			  "rawAngles" => 0,
                          "numCenter" => 0,
			  "numCluster" => 0,
			  "decisions" => 0,
			  "nonModels" => 0
                          );

sub new
  {
  my $class = shift @_;
  my $self = { @defaultDataMembers, @_ };
  bless $self, ref $class || $class;

  $self->readPDB($self->{element}) if ($self->{pathsFile} ne "" && $self->{element} ne "");
  return $self;
  }


# read PDB entries to create shells
sub readPDB
  {
  my $self = shift @_;
  my $element = shift @_;

  my $allShells = []; 
  open (PATH, $self->{pathsFile}); 
  foreach my $file (<PATH>)
    {
    chomp $file;

    $file =~ s/[\r\n]+$//;
    my $pdb = PDBEntry->new("singlePdbFile" => $file);
    my $atoms = $pdb->{atoms};

    my $shellsOfOnePDB = AtomShell->createShells($element, $atoms);
   
    ## Calculating number of zinc clusters
    my $zn = scalar (grep {$_->{"element"} eq $element} (@$atoms)); 
    my $cluster = $zn - scalar @$shellsOfOnePDB;
    $self->{numCenter} += $zn;
    $self->{numCluster} += $cluster;
    # print $$atoms[0]->{PDBid}, ": $zn zincs\n" if $zn > 15;

    my $residues = $pdb->{residues};
    my $sequences = $pdb->{sequences};
    foreach my $oneShell (@$shellsOfOnePDB) 
      {
      ## if an atom in the shell is not standard aa, mark its closest aa ligand
      my $closestAA = {};
      foreach my $lig (@{$oneShell->{shell}})
	{
	if (&Sequence::_aaCode($lig->{residueName}))
	  { $$closestAA{$lig} = $lig; }
	else 
	  { $$closestAA{$lig} = $lig->closest($atoms); }
	}
      $oneShell->{closestAA} = $closestAA;

      ## prepare residues for later analysis
      foreach my $ligand (@{$oneShell->{shellsOfOnePDB}})
	{
	map {$ligand->{chiAngle} = $_->chiOneAngle();} (grep {$ligand->resID() eq $_->residueID();} (@$residues));
	}

      ## store all chain-sequence pair of within its PDB into each shell
      $oneShell->{seqsOfPDB} = $sequences; #if (! $oneShell->{sequence});
      }

    push @$allShells, @$shellsOfOnePDB;
    }

  $self->{shells} = $allShells;
  }


## Bootstapping process for initializing coordination classification
sub bootstrapCoordination
  {
  my $self = shift @_;
  my $statOutFileName = (@_)? shift @_: "stats";

  my $stats = {};
  $self->calcDeviationCoordination();  ## 

  $$stats{"angle"} = $self->calcAngleStats();
  $$stats{"distance"} = $self->calcDistStats();
  $self->{stats} = $stats;

  &writeTableFile("$statOutFileName.0.txt", $stats);
  $self->printStats();
  }


## Expectation-Maximazation process for finer coordination classification
sub IAcoordination
  {
  my $self = shift @_;
  my $statOutFileName = shift @_;
  my $control = shift @_;
  my $threshold = shift @_;
  my $currStats = (@_)? (shift @_) : ($self->{stats});

  my $i = 1;
  my $oldStats = {};
  while ( $self->compareStats($oldStats) )
    {
    $oldStats = clone($currStats);
    $self->{coordinations} = $self->calcChiCoordination($control, $threshold);
    if ($control eq "n") {$self->calNonModel($threshold);}

    $currStats = {};
    $$currStats{"angle"} = $self->calcAngleStats();
    $$currStats{"distance"} = $self->calcDistStats();

    $self->{stats} = $currStats;
    $self->printStats($i);

    &writeTableFile("$statOutFileName.$i.txt", $currStats);

    if ($i == 50)
      { 
      print "Failed to stabilize\n";
      last; 
      }

    $i++;
    }
  }


##
sub calNonModel
  {
  my $self = shift @_;
  my $threshold = shift @_;

  my $stats;
  $$stats{"angle"} = $self->calcAngleStats();
  $$stats{"distance"} = $self->calcDistStats();

  my $models = $self->{coordinations};
  my $coordinations = {};
  my $shells = [];
  my $nonModels = []; 

  foreach my $model (keys %$models)
    {
    my $coord = $$models{$model};

    foreach my $struct (@$coord)
      {
      ## Seems to be error, leave it here for now.
      #$struct->bestTestStatistic("non", $stats);
    
      if ($struct->{bestCombo}->{probability} > $threshold)
        {
	push @{$$coordinations{$model}}, $struct;
        push @$shells, $struct->{shellObj}
	}
      else 
	{
	$struct->{class} = ref $struct;
	$struct->{angleList} = [$struct->angleList()];
	push @$nonModels, $struct;
	}
      }
    }

  if ($self->{nonModels})
    { push @{$self->{nonModels}}, @$nonModels; }
  else 
    { $self->{nonModels} = $nonModels;}
  $self->{coordinations} = $coordinations;
  $self->{shells} = $shells;

  }


## Use distance-only statistics in the chi probability test to determine how many ligands should a center take.
sub bindShellViaDist
  {
  my $self = shift @_;
  my $statOutFileName = (@_)? shift @_: "stats";
  my $stats = (@_)? (shift @_) : ($self->{stats});

  my ($tetree, $four, $five, $six, $more);
  foreach my $shell (@{$self->{shells}})
    {
      my $tpl = TrigonalPlanar->new(shellObj => $shell);
      my $tet = Tetrahedral->new(shellObj => $shell);
      my $tbp = TrigonalBipyramidal->new(shellObj => $shell);
      my $oct = Octahedral->new(shellObj => $shell);

      $tpl->bestDistChi($stats);
      $tet->bestDistChi($stats);
      $tbp->bestDistChi($stats);
      $oct->bestDistChi($stats);

    if (! defined $tpl->{bestCombo}) ## less than 3 ligands
      {
      next;
      }
    elsif (! defined $tet->{bestCombo}) ## 3 ligands
      {
      push @$tetree, $tpl;
      }
    elsif (! defined $tbp->{bestCombo}) ## 4 ligands
      {
      push @$four, $tet;
      }
    elsif (! defined $oct->{bestCombo}) ## 5 ligands
      {
      my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($tet, $tbp)))[0];
      if (ref $bestModel eq "Tetrahedral")
        { push @$four, $tet;}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { push @$five, $tbp;}
      }
    else
      {
      my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($tet, $tbp, $oct)))[0];
      if (ref $bestModel eq "Tetrahedral")
        { push @$four, $tet;}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { push @$five, $tbp;}
      elsif (ref $bestModel eq "Octahedral")
        { push @$six, $oct;}
      }
 
    ## print chi probabilities
    #print $shell->znID(), ": ";
    #print "th, ", $tet->{bestCombo}->{probability}, "; ";
    #print "tb, ", $tbp->{bestCombo}->{probability}, "; ";
    #print "oh, ", $oct->{bestCombo}->{probability}, "; ";
    #print "\n";
    }

  &writeTableFile("$statOutFileName.dist.txt", $stats);

  my $bestDist = {};
  %$bestDist = (  "three" => $tetree,
		  "four" => $four,
		  "five" => $five,
		  "six" => $six );

  $self->{coordinations} = $bestDist ;
  }


## After bestDistance, print out sequences in fasta format
sub printSequences
  {
  my $self = shift @_;
  my $outFile = shift @_;
  my $seqType = shift @_;
  my $headerType = shift @_;
  my $ligNum = shift @_;

  open (my $fileH, ">", $outFile) or die $!;
  foreach my $model (@{$self->{coordinations}{$ligNum}}) 
    {
    my $seqsOfPDB = $model->{shellObj}->{seqsOfPDB};
    my $znId = $model->{shellObj}->znID();

    my @headerLigs;
    if ($headerType eq "b")  ## original ligands
      { @headerLigs = @{$model->{bestCombo}->{ligands}}; }
    elsif ($headerType eq "s") ## shell ligands
      { @headerLigs = @{$model->{shellObj}->{shell}}; }
    elsif ($headerType eq "ss")  ### sub non-aa with aa in second shell but not binding ligands
      {
      my $secShell = [grep {my $temp = $_; ! grep {$_ eq $temp} (@{$model->{bestCombo}->{ligands}}); } (@{$model->{shellObj}->{secondShell}})];
      foreach my $lig (@{$model->{bestCombo}->{ligands}})
	{
	if (&Sequence::_aaCode($lig->{residueName}))
	  {push @headerLigs, $lig;}
	else
	  { 
	  if ($lig->closest($secShell)) 
	    {push @headerLigs, $lig->closest($secShell); }
          else
	    {print "flag, ", $lig->resID();} 
	  }
	} 
      }
    elsif ($headerType eq "c")  ### replace non-aa by closest aa
      {
      my $closestAA = $model->{shellObj}->{closestAA}; 
      @headerLigs = (map { $$closestAA{$_}; } (@{$model->{bestCombo}->{ligands}}));
      }

    my @ligId = (map {$_->resID().".".$_->{residueName}.".".$_->{atomName}.".".$_->{element} ;} (@headerLigs)); 
    
    my %chains;
    foreach my $lig (@headerLigs)
      {
      next if (! &Sequence::_aaCode($lig->{residueName}) ); #non-aa ligands do not count
      $chains{$lig->{chainID}} = 1;
      }

    if (scalar (keys %chains) == 0)
        {print "No protein ligands, $znId!\n";} 

    my $seqOfChain;
    my $seqs; 
    if ($seqType eq "s")
      {$seqs = $seqsOfPDB->{seqres};}
    elsif ($seqType eq "a")
      {$seqs = $seqsOfPDB->{atom};}
    elsif ($seqType eq "n")
      {$seqs = $seqsOfPDB->{number};}
   
    foreach my $ch (keys %chains)
      {
      foreach my $seq (@$seqs)
        {
        if ( substr((split('\|', $seq->{header}))[-1], 0, 1) eq $ch ) 
          { $seqOfChain = $seq; last; }
        }

      $seqOfChain->updateHeader($znId, @ligId);
      $seqOfChain->printFasta($fileH);
      }		
    }  

  close $fileH
  }


## compare statistics
sub compareStats
  {
  my $self = shift @_;
  my $oldStats = shift @_;
  my $cutoff = (@_)? shift @_: 0.0001;

  return 1 if (! keys %$oldStats);
  my $currStats = $self->{stats};

  foreach my $distOrAng (keys %$currStats)
    {
    return 2 if (! exists $$oldStats{$distOrAng} );

    foreach my $coordination (keys %{$$currStats{$distOrAng}})
      {
      if ($coordination eq "variance" || $coordination eq "standardDeviation")
        { return 3 if (abs ($$currStats{$distOrAng}{$coordination} - $$oldStats{$distOrAng}{$coordination}) > $cutoff) ; next;}

      return 4 if (! exists $$oldStats{$distOrAng}{$coordination} ); 

      my $curr = $$currStats{$distOrAng}{$coordination};
      my $old = $$oldStats{$distOrAng}{$coordination};
      if ($distOrAng eq "distance")
	{
	return 5 if ( grep { abs ($$curr{$_} - $$old{$_}) > $cutoff;} (keys %$curr) );
        return 6 if ( grep { abs ($$curr{$_} - $$old{$_}) > $cutoff;} (keys %$old) );
	}
      else
	{
        foreach my $angle (keys %$curr)
    	  {
	  return 7 if (! exists $$old{$angle});

	  return 8 if ( grep { abs ($$curr{$angle}{$_} - $$old{$angle}{$_}) > $cutoff; } (keys %{$$curr{$angle}}) );
          return 9 if ( grep { abs ($$curr{$angle}{$_} - $$old{$angle}{$_}) > $cutoff; } (keys %{$$old{$angle}}) );
	  }
	}
      }
    }

  foreach my $distOrAng (keys %$oldStats)
    {
    return 11 if (! exists $$currStats{$distOrAng} );

    foreach my $coordination (keys %{$$oldStats{$distOrAng}})
      {
      if ($coordination eq "variance" || $coordination eq "standardDeviation")
        { return 12 if ( abs ($$currStats{$distOrAng}{$coordination} - $$oldStats{$distOrAng}{$coordination}) > $cutoff ) ; next;}

      return 13 if (! exists $$currStats{$distOrAng}{$coordination} );

      my $curr = $$currStats{$distOrAng}{$coordination};
      my $old = $$oldStats{$distOrAng}{$coordination};
      if ($distOrAng eq "distance")
        {
        return 14 if ( grep { abs ($$curr{$_} - $$old{$_}) > $cutoff;} (keys %$curr) );
        return 15 if ( grep { abs ($$curr{$_} - $$old{$_}) > $cutoff;} (keys %$old) );
        }
      else
        {
        foreach my $angle (keys %$old)
          {
          return 16 if (! exists $$curr{$angle});

          return 17 if ( grep { abs ($$curr{$angle}{$_} - $$old{$angle}{$_}) > $cutoff; } (keys %{$$curr{$angle}}) );
          return 18 if ( grep { abs ($$curr{$angle}{$_} - $$old{$angle}{$_}) > $cutoff; } (keys %{$$old{$angle}}) );
          }
        }
      }
    }

  return 0;
  }

our $cgRelations = [
  {"name" => "Tetrahedral", "num" => 4, "parents" => [], "children" => ["TetrahedralV"], "siblings" => []},
  {"name" => "TrigonalBipyramidal", "num" => 5, "parents" => [], "children" => ["TrigonalBipyramidalVA", "TrigonalBipyramidalVP", "TrigonalPlanar"], , "siblings" => []},
  {"name" => "Octahedral", "num" => 6, "parents" => [], "children" => ["SquarePyramidalV", "SquarePlanar", "SquarePyramidal"], "siblings" => []},
  {"name" => "TrigonalPlanar", "num" => 3, "parents" => ["TrigonalBipyramidalVA", "TrigonalBipyramidal"], "children" => [], "siblings" => []},
  {"name" => "TetrahedralV", "num" => 3, "parents" => ["Tetrahedral"], "children" => [], "siblings" => []},
  {"name" => "TrigonalBipyramidalVA", "num" => 4, "parents" => ["TrigonalBipyramidal"], "children" => ["TrigonalPlanar"], "siblings" => ["TrigonalBipyramidalVP"]},
  {"name" => "TrigonalBipyramidalVP", "num" => 4, "parents" => ["TrigonalBipyramidal"], "children" => [], "siblings" => ["TrigonalBipyramidalVA"]},
  {"name" => "SquarePyramidalV", "num" => 4, "parents" => ["SquarePyramidal", "Octahedral"], "children" => [], "siblings" => ["SquarePlanar"]},
  {"name" => "SquarePlanar", "num" => 4, "parents" => ["SquarePyramidal", "Octahedral"], "children" => [], "siblings" => ["SquarePyramidalV"]},
  {"name" => "SquarePyramidal", "num" => 5, "parents" => ["Octahedral"], "children" => ["SquarePlanar", "SquarePyramidalV"], "siblings" => []}
]; 

## Classify coordination using chi statistics
sub calcChiCoordination
  {
  my $self = shift @_;
  my $control = shift @_;
  my $threshold = shift @_;
  my $stats = (@_)? (shift @_) : ($self->{stats});

  my $decisions = {};
  my $coordinations = {};
  foreach my $shell (@{$self->{shells}})
    {
    # major coordinations
    my $tet = Tetrahedral->new(shellObj => $shell);
    my $tbp = TrigonalBipyramidal->new(shellObj => $shell);
    my $oct = Octahedral->new(shellObj => $shell);

    # 3 ligands
    my $tpl = TrigonalPlanar->new(shellObj => $shell);
    my $tev = TetrahedralV->new(shellObj => $shell);

    # 4 ligands, trignal bipyramidal related
    my $bva = TrigonalBipyramidalVA->new(shellObj => $shell);
    my $bvp = TrigonalBipyramidalVP->new(shellObj => $shell);

    # 4 ligands, square pyramidal related
    my $spv = SquarePyramidalV->new(shellObj => $shell);
    my $spl = SquarePlanar->new(shellObj => $shell);

    # 5 ligands
    my $spy = SquarePyramidal->new(shellObj => $shell);

    $tet->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $tbp->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $oct->bestTestStatistic("chi", $control, $threshold, 0, $stats);

    $tpl->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $tev->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $bva->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $bvp->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spv->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spl->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spy->bestTestStatistic("chi", $control, $threshold, 0, $stats);

    my @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} ($tet, $tbp, $oct, $bva, $bvp, $spv, $spl, $spy)));

#print "\n", $shell->znID(), "\n";
#map {print ref $_, ": ", $_->{bestCombo}->{probability}, "\n"} (@models);

    my $maxNum;
    my $unusables;
    foreach my $model (@models)
      {
      my $relation = (grep {$$_{"name"} eq ref $model} (@$cgRelations))[0];
      my $num = $$relation{"num"};
      if ($num > $maxNum) {$maxNum = $num;}
      }

    if ($maxNum < 4) ## less than 4 ligands
      {
      if (! defined $tev->{bestCombo} && ! defined $tpl->{bestCombo}) 
	{ $$decisions{"012"}++; }
      else 
        {
        my @mods = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} ($tev, $tpl)));
        if ($mods[0]->{bestCombo}->{probability} > 0.01)
          {
          my $modRef = ref $mods[0];
          push @{$$coordinations{$modRef}}, $mods[0];
          $$decisions{"3.". $modRef} += 1;
          }
        else
          {$$decisions{"3.None"} += 1;}
        next;
        }
      }
    else ## 4, 5, 6 ligands
      {
      ## set the probability threshold, remove low prob ones from statistics calculation. 
      if ($control eq "p" && $models[0]->{bestCombo}->{probability} < $threshold) 
	{
	push @$unusables, $models[0];
        $$decisions{$maxNum. ".Unusable"} += 1;
	next;
	}

      if (ref $models[0] eq "TrigonalBipyramidalVA" && ref $models[1] eq "Tetrahedral" && ($models[0]->{bestCombo}->{probability} < (2 * $models[1]->{bestCombo}->{probability})))
        {
        my $modelRef = ref $models[1];
	push @{$$coordinations{$modelRef}}, $models[1];
	$$decisions{$maxNum. ".".  $modelRef} += 1;
#print "model class: $modelRef\n";
	}
      else
	{
	my $dec = $maxNum;
	my @track;
	for (my $i = 0; $i <= $#models; $i++)
	  {
          my $modelRef = ref $models[$i];
	  my $relation = (grep {$$_{"name"} eq $modelRef} (@$cgRelations))[0];

	  if ( grep {ref $models[$i+1] eq $_;} (@{$$relation{"parents"}}) )
	    { push @track, "p"; }
          elsif ( grep {ref $models[$i+1] eq $_;} (@{$$relation{"siblings"}}, @{$$relation{"children"}}) )
            { push @track, "s"; }
	  else
	    {
#print "track: ", @track, "\n";
	    if ($i == 0)
	      {	
              push @{$$coordinations{$modelRef}}, $models[$i];
	      $dec = $dec. ".". $modelRef;
              $$decisions{$dec} += 1;
	      }
	    else
	      {
	      while ($track[-1] eq "s") {pop @track;}
	      my $maxInd = @track;

	      $modelRef = ref $models[$maxInd];
              push @{$$coordinations{$modelRef}}, $models[$maxInd];

	      map {$dec = $dec.".".ref $models[$_]} (0..$maxInd) ;
              $$decisions{$dec} += 1;
	      }
#print "track: ", @track, "\n";
#print "decision: $dec\n";
#print "model class: $modelRef\n";
	    last;
	    }
	  }
	}
      }
    }

  $self->{decisions} = $decisions ;
  $self->{coordinations} = $coordinations ;
  }


## Classify coordination using deviation
sub calcDeviationCoordination
  {
  my $self = shift @_;

  my $shells = $self->{shells};
  my ($tetrahedrals, $trigonalbipyramidals, $octahedrals);
  foreach my $shell (@$shells)
    {
    my $tet = Tetrahedral->new(shellObj => $shell);
    my $tbp = TrigonalBipyramidal->new(shellObj => $shell);
    my $oct = Octahedral->new(shellObj => $shell);

    $tet->bestTestStatistic("dev");
    $tbp->bestTestStatistic("dev");
    $oct->bestTestStatistic("dev");

    my $bestModel = (sort {$a->{bestCombo}->{deviation} <=> $b->{bestCombo}->{deviation}} (grep {$_->{bestCombo}->{deviation} != 0;} ($tet, $tbp, $oct)))[0];
    if (ref $bestModel eq "Tetrahedral")
      { push @$tetrahedrals, $tet; }
    elsif (ref $bestModel eq "TrigonalBipyramidal")
      { push @$trigonalbipyramidals, $tbp; }
    elsif (ref $bestModel eq "Octahedral")
      { push @$octahedrals, $oct; }

    }

  my $coordinations = {};
  %$coordinations = ( "tetrahedrals" => $tetrahedrals,
                     "trigonalBipyramidals" => $trigonalbipyramidals,
                     "octahedrals" => $octahedrals );

  $self->{coordinations} = $coordinations;

  }


sub printStats
  {
  my $self = shift @_;
  my $i = (@_)? shift @_ : 0;

  if ($i ==0)
    { print "Deviation sorted:" ;}
  else
    { print "Chi sorted, i = $i:" ;}

  foreach my $coord (sort keys %{$self->{coordinations}})
    {
    print "\n$coord = ", scalar @{$self->{coordinations}->{$coord}} if $self->{coordinations}->{$coord};
    }
  print "\n\n";
  }


## Given coordination classification, calculated angle statistics
sub calcAngleStats
  {
  my $self = shift @_;

  my $coordSets = $self->{coordinations};
  my $angleStats = {};

  my $coordinationAgnles = {};
  foreach my $coord (keys %$coordSets)
    {
    my $coordSet = $$coordSets{$coord};
    next if (! $coordSet);
    map { $_->calcExpectedAngleStats($coordinationAgnles) ;} (@$coordSet);
    }

  my $totalDev = 0;
  my $totaln = 0;
  foreach my $coordination (keys %$coordinationAgnles)
    {
    foreach my $angle (keys %{$$coordinationAgnles{$coordination}})
      {
      my $stats = RawStatistics->new("variables" => $$coordinationAgnles{$coordination}{$angle} ) ;
      my $mean = $stats->calcMean() ;
      my $deviation = $stats->calcDeviation() ;
      my $count = $stats->count();

      if ($coordination eq "squarePlanar" || $angle ne "180")
	{
        $totalDev += $deviation;
        $totaln += $count;
	}

      $$angleStats{$coordination}{$angle}{"mean"} = $mean;
      $$angleStats{$coordination}{$angle}{"count"} = $count;
      #$$angleStats{$coordination}{$angle}{"min"} = $stats->min();
      #$$angleStats{$coordination}{$angle}{"max"} = $stats->max();
      $$angleStats{$coordination}{$angle}{"variance"} = $stats->calcVariance;
      $$angleStats{$coordination}{$angle}{"standardDeviation"} = $stats->calcStd();
      }
    }
 
  my $variance = $totalDev/$totaln ;
  $$angleStats{"variance"} = $variance;
  $$angleStats{"standardDeviation"} = $variance ** 0.5;

  $self->{rawAngles} = $coordinationAgnles;
  return $angleStats;
  }

## Given coordination classification, calculated distance statistics       
sub calcDistStats
  {
  my $self = shift @_;

  my $coordSets = $self->{coordinations};
  my $elementDists = {};
  my $distStats = {};

  foreach my $coord (keys %$coordSets)
    {
    my $coordSet = $$coordSets{$coord};

    foreach my $model (@$coordSet)
      {
      foreach my $ligand ( @{$model->{bestCombo}->{ligands}} )
        {
        my $element = $ligand->{element};
        push (@{$$elementDists{$element}}, $model->{shellObj}->{center}->distance($ligand));
        }
      }
    }

  foreach my $element (keys %$elementDists)
    {
    map {push (@{$$elementDists{"average"}}, $_) ;} (@{$$elementDists{$element}});

    if (@{$$elementDists{$element}} > 30)
      {
      my $stats = RawStatistics->new("variables" => $$elementDists{$element} ) ;
      my $mean = $stats->calcMean() ;
      my $var =+ $stats->calcVariance() ;

      $$distStats{$element} = {"mean" => $mean, "variance" => $var, "count" => $stats->count(), "max" => $stats->max(), "min" => $stats->min(), "standardDeviation" => $stats->calcStd()};
      }
    }

  my $stats = RawStatistics->new("variables" => $$elementDists{"average"} ) ;
  my $mean = $stats->calcMean() ;
  my $var =+ $stats->calcVariance() ;

  $$distStats{"average"} = {"mean" => $mean, "variance" => $var, "count" => $stats->count(), "max" => $stats->max(), "min" => $stats->min(), "standardDeviation" => $stats->calcStd()};

  return $distStats;
  }




### Obsolete, see methods in RawStatistics
sub calcMeanVarMaxMin
  {
  if (@_)
    {
    my $mean = 0;
    my $var = 0;
    my $max = $_[0];
    my $min = $_[0];

    foreach my $value (@_)
      {
      $mean += $value;
      if ($value < $min)
        { $min = $value; }
      elsif ($value > $max)
        { $max = $value; }
      }

    $mean = $mean / scalar(@_);

    foreach my $value (@_)
      { $var += ($mean - $value) ** 2; }

    $var = $var / scalar(@_);

    return (scalar(@_),$mean,$var,$max,$min);
    }

  return (0,0,0,0,0)
  }












 
