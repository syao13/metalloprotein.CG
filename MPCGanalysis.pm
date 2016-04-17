##!/usr/bin/perl 
#===============================================================================
#
#         FILE:  MPCGanalysis.pm
#
#        USAGE:  ./MPCGanalysis.pm 
#
#  DESCRIPTION:  Metalloprotein (MP) coordination geometry (CG) analysis,
#  		 A comprehensive process object that makes uses of most of the
#  		 other objects and conducts the metal coordination classification 
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

package MPCGanalysis;
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
use PentagonalBipyramidal;
use Cube;
use SquareAntiprismaticMonocapped;
use SquareAntiprismaticBicapped;
use PentagonalBipyramidalVA;
use PentagonalBipyramidalVP;
use SquareAntiprismatic;
use SquareAntiprismaticV;
use TrigonalPrismatic;
use TrigonalPrismaticV;
use HexagonalBipyramidal;
use HexagonalBipyramidalVA;
use HexagonalBipyramidalVP;

use threads;
use Thread::Queue;
use Clone 'clone';
#use Data::Dumper::Concise; # human readable, code in iaCoordination
use JSON; # For other programs to read back, code in iaCoordination 
use JSON -convert_blessed_universally;
use Time::HiRes qw(time);
use POSIX qw(strftime);

our @defaultDataMembers = (
                          "pathsFile" => 0,  # a file that contains pdb-file paths
			  "element" => 0,
			  "majorCGs" => 0,
			  "minLigNum" => 0,
			  "shells" => 0,  # ref to AtomShell object
			  "coordinations" => 0,  
			  "stats" => 0,
			  "rawAngles" => 0,
                          "numCenter" => 0,
			  "numCluster" => 0,
			  "decisions" => 0,
			  "nonModels" => 0,
			  "unusable" => 0,
			  "usable" => 0,
			  "shellOpt" => 0,
                          "shellCutoff" => 0
                          );

our $cgRelations = [
  {"name" => "Tetrahedral", "num" => 4, "parents" => [], "children" => ["TetrahedralV"], "siblings" => []},
  #{"name" => "TetrahedralV", "num" => 3, "parents" => ["Tetrahedral"], "children" => [], "siblings" => []},

  {"name" => "TrigonalBipyramidal", "num" => 5, "parents" => [], "children" => ["TrigonalBipyramidalVA", "TrigonalBipyramidalVP", "TrigonalPlanar"], "siblings" => []},
  {"name" => "TrigonalBipyramidalVA", "num" => 4, "parents" => ["TrigonalBipyramidal"], "children" => ["TrigonalPlanar"], "siblings" => ["TrigonalBipyramidalVP"]},
  {"name" => "TrigonalBipyramidalVP", "num" => 4, "parents" => ["TrigonalBipyramidal"], "children" => [], "siblings" => ["TrigonalBipyramidalVA"]},
  #{"name" => "TrigonalPlanar", "num" => 3, "parents" => ["TrigonalBipyramidalVA", "TrigonalBipyramidal"], "children" => [], "siblings" => []},

  {"name" => "Octahedral", "num" => 6, "parents" => [], "children" => ["SquarePyramidalV", "SquarePlanar", "SquarePyramidal"], "siblings" => []},
  {"name" => "SquarePyramidal", "num" => 5, "parents" => ["Octahedral"], "children" => ["SquarePlanar", "SquarePyramidalV"], "siblings" => []},
  {"name" => "SquarePyramidalV", "num" => 4, "parents" => ["SquarePyramidal", "Octahedral"], "children" => [], "siblings" => ["SquarePlanar"]},
  {"name" => "SquarePlanar", "num" => 4, "parents" => ["SquarePyramidal", "Octahedral"], "children" => [], "siblings" => ["SquarePyramidalV"]},

  {"name" => "TrigonalPrismatic", "num" => 6, "parents" => [], "children" => ["TrigonalPrismaticV"], "siblings" => []},
  {"name" => "TrigonalPrismaticV", "num" => 5, "parents" => ["TrigonalPrismatic"], "children" => [], "siblings" => []},

  {"name" => "PentagonalBipyramidal", "num" => 7, "parents" => [], "children" => ["PentagonalBipyramidalVA", "PentagonalBipyramidalVP"], "siblings" => []},
  {"name" => "PentagonalBipyramidalVA", "num" => 6, "parents" => ["PentagonalBipyramidal"], "children" => [], "siblings" => ["PentagonalBipyramidalVP"]},
  {"name" => "PentagonalBipyramidalVP", "num" => 6, "parents" => ["PentagonalBipyramidal"], "children" => [], "siblings" => ["PentagonalBipyramidalVA"]},

  {"name" => "SquareAntiprismatic", "num" => 8, "parents" => [], "children" => ["SquareAntiprismaticV"], "siblings" => []},
  {"name" => "SquareAntiprismaticV", "num" => 7, "parents" => ["SquareAntiprismatic"], "children" => [], "siblings" => []},

  {"name" => "HexagonalBipyramidal", "num" => 8, "parents" => [], "children" => ["HexagonalBipyramidalVA", "HexagonalBipyramidalVP"], "siblings" => []},
  {"name" => "HexagonalBipyramidalVA", "num" => 7, "parents" => ["HexagonalBipyramidal"], "children" => [], "siblings" => ["HexagonalBipyramidalVP"]},
  {"name" => "HexagonalBipyramidalVP", "num" => 7, "parents" => ["HexagonalBipyramidal"], "children" => [], "siblings" => ["HexagonalBipyramidalVA"]},

  {"name" => "SquareAntiprismaticMonocapped", "num" => 9, "parents" => [], "children" => [], "siblings" => []},
  {"name" => "SquareAntiprismaticBicapped", "num" => 10, "parents" => [], "children" => [], "siblings" => []}
]; 


## The combined lm slope of bond length std vs. resolution
our $slope = 0.0562;

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

#my $now_string = localtime;
#print "start, $now_string", "\n";
  while (my $file = <PATH>)
    {
    chomp $file;
    my $pdb = PDBEntry->new("singlePdbFile" => $file, "metal" => $element);
    my $atoms = $pdb->{atoms};
#my $now_string = localtime;
#print "\n", substr($file,42,4), ",", $pdb->{symNum}, ", $now_string\n";
    my $shellsOfOnePDB = ($self->{shellCutoff})? AtomShell->createShells($element, $atoms, 1.3, $self->{shellCutoff}, $self->{shellElement}) : AtomShell->createShells($element, $atoms);

#my $now_string = localtime;
#print "after getting shells, $now_string\n";
   
    ## Calculating number of zinc clusters
    my $metal = scalar (grep {$_->{"element"} eq $element} (@$atoms)); 
    my $cluster = $metal - scalar @$shellsOfOnePDB;
    $self->{numCenter} += $metal;
    $self->{numCluster} += $cluster;
    # print $$atoms[0]->{PDBid}, ": $metal zincs\n" if $metal > 15;

    my $residues = $pdb->{residues};
    my $sequences = $pdb->{sequences};
    foreach my $oneShell (@$shellsOfOnePDB) 
      {
      if (scalar @{$oneShell->{shell}} > 3) { $self->{usable} += 1; }
      else { $self->{unusable} += 1; }

      ## if an atom in the shell is not standard aa, mark its closest aa ligand
      ## commented out due to taking too long to run. 16.4.15
      #my $closestAA = {};
      #foreach my $lig (@{$oneShell->{shell}})
	#{
	#if (&Sequence::_aaCode($lig->{residueName}))
	#  { $$closestAA{$lig} = $lig; }
	#else 
	#  { $$closestAA{$lig} = $lig->closest($atoms); }
	#}
      #$oneShell->{closestAA} = $closestAA;

      ## prepare residues for later analysis
      #foreach my $ligand (@{$oneShell->{shellsOfOnePDB}})
	#{
	#map {$ligand->{chiAngle} = $_->chiOneAngle();} (grep {$ligand->resID() eq $_->residueID();} (@$residues));
	#}

      ## store all chain-sequence pair of within its PDB into each shell
      $oneShell->{seqsOfPDB} = $sequences; #if (! $oneShell->{sequence});
      }

    push @$allShells, @$shellsOfOnePDB;
    }

#my $now_string = localtime;
#print "end, $now_string", "\n";

  close PATH;
  $self->{shells} = $allShells;
  }


## Bootstapping process for initializing coordination classification
sub bootstrapCoordination
  {
  my $self = shift @_;
  my $statOutFileName = (@_)? shift @_: "stats";
  $statOutFileName = "stats" unless ($statOutFileName);

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
    $oldStats = clone($currStats); ## a nested copy
    $self->calcChiCoordination($control, $threshold);
    if ($control eq "n") {$self->calNonModel($threshold);}

    $currStats = {};
    $$currStats{"angle"} = $self->calcAngleStats();
    $$currStats{"distance"} = $self->calcDistStats();

    $self->{stats} = $currStats;
    $self->printStats($i);

    &writeTableFile("$statOutFileName.$i.txt", $currStats);

    if ($i == 20)
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

  my %numToLet = ( 2 => "two", 3 => "three", 4 => "four", 5 => "five", 6 => "six", 7 => "seven", 8 => "eight", 9 => "nine", 10 => "ten");
  my $coordinations = {};
  foreach my $shell (@{$self->{shells}}) 
    {
    my @models;
    #foreach my $cg ("TrigonalPlanar", @{$self->{majorCGs}}, "PentagonalBipyramidal", "Cube", "SquareAntiprismaticMonocapped", "SquareAntiprismaticBicapped")
    foreach my $cg (@{$self->{majorCGs}}, "SquareAntiprismatic", "SquareAntiprismaticMonocapped", "SquareAntiprismaticBicapped")
      {
      my $cgObj = $cg->new("shellObj" => $shell);
      $cgObj->bestDistChi($stats);
      push @models, $cgObj;
      }

    my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} (@models)))[0];
    next if (! $bestModel);

    my $relation = (grep {$$_{"name"} eq ref $bestModel} (@$cgRelations))[0];
    my $numLig = $numToLet{$$relation{"num"}};
    push @{$$coordinations{$numLig}}, $bestModel;

    # print chi probabilities
    #print $shell->metalID(), "; ";
    #map {print $_->{bestCombo}->{probability}, "; "} (@models);
    #print "\n";
    #print "$bestModel\n";
    #print "\n";
    }

  $self->{coordinations} = $coordinations;
  }

sub shellViaAdjustDistStd
  {
  my $self = shift @_;
  my $statOutFileName = (@_)? shift @_: "stats";
  my $stats = (@_)? (shift @_) : ($self->{stats});
  my $blStats = $$stats{"distance"};

  my $coordinations = {};
  foreach my $shell (@{$self->{shells}})
    {
    my $finalShell = [];

    ## Get the best alternate location.
    my $alternate = {};
    my @altAtoms;
    foreach my $ligand (@{$shell->{shell}})
      { 
      push (@altAtoms, $ligand->resID()) if ($$alternate{$ligand->resID()} && $$alternate{$ligand->resID()}{$ligand->{alternateLocation}} != 1);
      $$alternate{$ligand->resID()}{$ligand->{alternateLocation}} = 1;
      }

    if (@altAtoms) ## If there are two alternate atoms
      {
      my $altByRes = {};
      foreach my $ligand (@{$shell->{shell}})
	{
	if (grep {(split(/\./, $_))[0] eq $ligand->resID()} (@altAtoms))
	  {
	  my $resolution = ($ligand->{resolution} == -1)? 2.5 : $ligand->{resolution};
          my $adjStd = ($resolution - $$blStats{$ligand->{element}}{resolutionAvg}) * $slope + $$blStats{$ligand->{element}}{standardDeviation};
          my $score = abs($shell->{center}->distance($ligand) - $$blStats{$ligand->{element}}{mean})/$adjStd;
	  $$altByRes{$ligand->resID()}{$ligand->{alternateLocation}}{"score"} += $score;
          push @{$$altByRes{$ligand->resID()}{$ligand->{alternateLocation}}{"shells"}}, $ligand;
	  }
	else
	  {
	  my $resolution = ($ligand->{resolution} == -1)? 2.5 : $ligand->{resolution};
          my $adjStd = ($resolution - $$blStats{$ligand->{element}}{resolutionAvg}) * $slope + $$blStats{$ligand->{element}}{standardDeviation};
          if ( abs($shell->{center}->distance($ligand) - $$blStats{$ligand->{element}}{mean}) <= $adjStd * 2.5 )
            { push @$finalShell, $ligand; }
	  }
 	}

      foreach my $res (keys %$altByRes) ## choose the best from all alternate locations based on each residue.
	{
	foreach my $altID (keys %{$$altByRes{$res}})
	  { $$altByRes{$res}{$altID}{"score"} = $$altByRes{$res}{$altID}{"score"} / (scalar @{$$altByRes{$res}{$altID}{"shells"}}); }
	my $minScore = (sort {$a <=> $b} (map {$$_{"score"}} (values %{$$altByRes{$res}})))[0];
	my $best = (grep {$$_{"score"} == $minScore} (values %{$$altByRes{$res}}))[0];
	push @$finalShell, @{$$best{"shells"}}
	}
      }	
    else ## No alternate atoms exist 
      {
      foreach my $ligand (@{$shell->{shell}})
        {
        my $resolution = ($ligand->{resolution} == -1)? 2.5 : $ligand->{resolution};
        my $adjStd = ($resolution - $$blStats{$ligand->{element}}{resolutionAvg}) * $slope + $$blStats{$ligand->{element}}{standardDeviation};
        if ( abs($shell->{center}->distance($ligand) - $$blStats{$ligand->{element}}{mean}) <= $adjStd * 2.5 )
	  { 
          push @$finalShell, $ligand;  
	  }
        }
      }

#print join(", ", "before", $shell->metalID(), map {$_->resID();} (@$finalShell)), "\n";
    my $excludeInd = 0; ## eliminate ligands with unreasonable atom-atom distances
    for (my $x=0; $x < @$finalShell -1; $x++)
      {
      last if ($excludeInd);
      for (my $y=$x+1; $y < @$finalShell; $y++)
        {
	last if ($excludeInd);
        my $distBtLigs = $$finalShell[$x]->distance($$finalShell[$y]) ;
        if ($distBtLigs < 1.5 || $distBtLigs > 6.0 )
	  { 
#print $$finalShell[$x]->resID(), ", ", $$finalShell[$x]->resID(), ", ", "$distBtLigs\n";
	  $excludeInd = 1 ; 
	  }
        }
      }

#print "tooclose, ", $shell->metalID(), ": ", join (", ", map { $_->atomID().".". $_->{alternateLocation} } (@$finalShell)), "\n" if ($excludeInd);
    next if ($excludeInd);
    foreach my $ligand (@$finalShell) ## eliminate symmetry related ligands that is not all water
      {
      if (substr ($ligand->{chainID}, 0, 1) eq "#" && $ligand->{residueName} ne "HOH")
	{ 
	$excludeInd = 1; 
	last;
	}
      }
#print "symmetry, ", $shell->metalID(), ": ", join (", ", map { $_->atomID().".". $_->{alternateLocation} } (@$finalShell)), "\n" if ($excludeInd);
    next if ($excludeInd);
    my @waters = grep {$_->{residueName} eq "HOH"} (@$finalShell); ## water cannot be the majority of the ligands
#print "water, ", $shell->metalID(), ": ", join (", ", map { $_->atomID().".". $_->{alternateLocation} } (@$finalShell)), "\n" if ((scalar @waters) * 2 > (scalar @$finalShell));
    next if ((scalar @waters) * 2 > (scalar @$finalShell));

#print join(", ", "after", $shell->metalID(), map {$_->resID();} (@$finalShell)), "\n";
    my $numLig = @$finalShell;
    next if ($numLig < 4 || $numLig > 10);

    my $cg =  (grep {$$_{"num"} eq $numLig;} (@$cgRelations))[0];
    my $shellObj = AtomShell->new("center" => $shell->{center}, "shell" => $finalShell, "seqsOfPDB" => $shell->{seqsOfPDB});
    my $cgObj = $$cg{"name"}->new("shellObj" => $shellObj);
    $cgObj->bestDistChi($stats);    

#print ref $cgObj, ", ";
#print $cgObj->{bestCombo}->{probability}, "\n\n"; 

    my %numToLet = ( 2 => "two", 3 => "three", 4 => "four", 5 => "five", 6 => "six", 7 => "seven", 8 => "eight", 9 => "nine", 10 => "ten");
    my $numLigL = $numToLet{$numLig};
    push @{$$coordinations{$numLigL}}, $cgObj;
    } 

  $self->{coordinations} = $coordinations;
  }


## After bestDistance, print out sequences in fasta format
sub printSequences
  {
  my $self = shift @_;
  my $outFile = shift @_;
  my $seqType = shift @_;
  my $headerType = shift @_;

  open (my $fileH, ">>", $outFile) or die $!;
  foreach my $ligNum (keys %{$self->{coordinations}})
    {
    foreach my $model (@{$self->{coordinations}{$ligNum}}) 
      {
      my $seqsOfPDB = $model->{shellObj}->{seqsOfPDB};
      my $metalId = $model->{shellObj}->metalID();

      my @headerLigs;
      if ($headerType eq "b")  ## original ligands
        { 
print $model->{shellObj}->metalID(), "\n" if (ref $model->{bestCombo}->{ligands} ne "ARRAY");
	@headerLigs = @{$model->{bestCombo}->{ligands}}; }
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

    #if (scalar (keys %chains) == 0)
    #    {print "No protein ligands, $metalId!\n";} 

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

        $seqOfChain->updateHeader($metalId, @ligId);
        $seqOfChain->printFasta($fileH);
        }		
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
      if ($coordination eq "variance" || $coordination eq "standardDeviation" || $coordination eq "pooledStd" || $coordination eq "pooledVar")
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
      if ($coordination eq "variance" || $coordination eq "standardDeviation" || $coordination eq "pooledStd" || $coordination eq "pooledVar")
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

## Classify coordination using chi statistics
sub calcChiCoordination
  {
  my $self = shift @_;
  my $control = shift @_;
  my $threshold = shift @_;
  my $stats = (@_)? (shift @_) : ($self->{stats});

  ## Acquire all CGs from major CGs
  #my %allCGs;
  #foreach my $major (@{$self->{majorCGs}})
  #  {
  #  $allCGs{$major} = 1;

  #  my $relation = (grep {$$_{"name"} eq $major } (@$cgRelations))[0];
  #  map {$allCGs{$_} = 1} (@{$$relation{"children"}});
  #  }
 
  ## all CGs in the cgRelations on top, regardless of the major ones passed in from bootstrap. 
  ## It is a array now compared to a hash as above.
  my @allCGs = map {$$_{"name"}} (grep {$$_{"num"} < 9 } (@$cgRelations));

  my $worker = sub 
    {
    my $tid = threads->tid;
    my ($shell, $Qwork, $Qresults ) = @_;
    while( my $work = $Qwork->dequeue() )
      {
      my $relation = (grep {$$_{"name"} eq $work } (@$cgRelations))[0];

print $work, ", ", $$relation{"num"}, ", ", $self->{minLigNum}, ", ", scalar @{$shell->{shell}}, "\n";
      next if ($$relation{"num"} < $self->{minLigNum} || $$relation{"num"} > @{$shell->{shell}});

      my $cgObj = $work->new("shellObj" => $shell);
      $cgObj->bestTestStatistic("chi", $control, $threshold, 0, $stats);
      #$Qresults->enqueue($work, $cgObj->{bestCombo}->{probability});
      $Qresults->enqueue($cgObj);

my $now_string = localtime;
print "$work, ", $cgObj->{bestCombo}->{probability}, ",  $now_string\n";
      }
    $Qresults->enqueue( undef ); ## Signal this thread is finished
    };

  my $decisions = {};
  my $coordinations = {};
  foreach my $ (@{$self->{shells}})
    {
    my $shell = $$self{"shells"}[$i];
print "\n", $shell->metalID(), "; ";
my $now_string = localtime;
print "$now_string\n";

    ## Parallel processing 
    my $THREADS = 8;
    my $Qwork = new Thread::Queue;
    my $Qresults = new Thread::Queue;

    foreach my $one (@allCGs)
      { $Qwork->enqueue($one); } #load the shared queue
    $Qwork->enqueue( (undef) x $THREADS ); # Tell the queue there are no more work items

    my @models;
    my @thread_pool = map { threads->create( $worker, $shell, $Qwork, $Qresults ) } (1 .. $THREADS);

my $leftW = $Qwork->pending();
my $leftR = $Qresults->pending();
print "before, $leftW, $leftR.\n";  

    while(my $result = $Qresults->dequeue())
      {push @models, $result;}
my $leftW = $Qwork->pending();
my $leftR = $Qresults->pending();
print "middle, $leftW, $leftR.\n";
    foreach my $th (@thread_pool)
      { $th->join; }
my $leftW = $Qwork->pending();
my $leftR = $Qresults->pending();
print "after, $leftW, $leftR.\n";

print join(", ", map {ref $_,  $_->{bestCombo}->{probability}} (@models)), ": before\n";
    @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} (@models)));
print join(", ", map {ref $_,  $_->{bestCombo}->{probability}} (@models)), ": after\n";

    ## Create all CG objects, single thread processing
#    my @models;
#    foreach my $cg (@allCGs)
#      {
#      my $relation = (grep {$$_{"name"} eq $cg } (@$cgRelations))[0];
#      next if ($$relation{"num"} < $self->{minLigNum} || $$relation{"num"} > @{$shell->{shell}});

#      my $cgObj = $cg->new(shellObj => $shell);
#      $cgObj->bestTestStatistic("chi", $control, $threshold, 0, $stats);
#      push @models, $cgObj;
    
#my $now_string = localtime;
#print "$cg, ", $cgObj->{bestCombo}->{probability}, ",  $now_string\n";
#      }
#    @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} (@models)));
    
    
    my $maxNum;
    my $unusables;
    ## Find the maximum number of ligands each metal structure has
    foreach my $model (@models)
      {
      my $relation = (grep {$$_{"name"} eq ref $model} (@$cgRelations))[0];
      my $num = $$relation{"num"};
      if ($num > $maxNum) {$maxNum = $num;}
      }

    if ($maxNum < $self->{minLigNum}) ## less than 4 ligands
      {
      my $tpl = TrigonalPlanar->new(shellObj => $shell);
      my $tev = TetrahedralV->new(shellObj => $shell);
      $tpl->bestTestStatistic("chi", $control, $threshold, 0, $stats);
      $tev->bestTestStatistic("chi", $control, $threshold, 0, $stats);

      if (! defined $tev->{bestCombo} && ! defined $tpl->{bestCombo}) 
	{ 
print "0;0;0\n";
	$$decisions{"012"}++; 
	}
      else 
        {
        my @mods = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{probability} != 0;} ($tev, $tpl)));
        if ($mods[0]->{bestCombo}->{probability} > 0.01)
          {
          my $modRef = ref $mods[0];
          push @{$$coordinations{$modRef}}, $mods[0];
          $$decisions{"3.". $modRef} += 1;
print "3lig.$modRef\n";
          }
        else
          {$$decisions{"3.None"} += 1;}
print $mods[0]->{bestCombo}->{probability}, "; ", $mods[1]->{bestCombo}->{probability}, "; 00\n";
        next;
        }
      }
    else ## 4, 5, 6 ligands
      {
print $models[0]->{bestCombo}->{probability}, "; 456\n";
      
      ## set the probability threshold, remove low prob ones from statistics calculation. 
      if ($control eq "p" && $models[0]->{bestCombo}->{probability} < $threshold) 
	{
	push @$unusables, $models[0];
        $$decisions{$maxNum. ".Unusable"} += 1;
	next;
	}

      ## Bva prob has to be more than twice of Tet to be considered Bva, otherwise Tet
      if (ref $models[0] eq "TrigonalBipyramidalVA" && ref $models[1] eq "Tetrahedral" && ($models[0]->{bestCombo}->{probability} < (2 * $models[1]->{bestCombo}->{probability})))
        {
        my $modelRef = ref $models[1];
	push @{$$coordinations{$modelRef}}, $models[1];
	$$decisions{$maxNum. ".".  $modelRef} += 1;

print $models[1]->{bestCombo}->{probability}, "; bva\n";
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
	    if ($i == 0)
	      {	
              push @{$$coordinations{$modelRef}}, $models[$i];
	      $dec = $dec. ".". $modelRef;
              $$decisions{$dec} += 1;
print $models[0]->{bestCombo}->{probability}, "; only $modelRef\n";
	      }
	    else
	      {
	      while ($track[-1] eq "s") {pop @track;}
	      my $maxInd = @track;

	      $modelRef = ref $models[$maxInd];
              push @{$$coordinations{$modelRef}}, $models[$maxInd];

	      map {$dec = $dec.".".ref $models[$_]} (0..$maxInd) ;
              $$decisions{$dec} += 1;
print $models[$maxInd]->{bestCombo}->{probability}, "; major CG to $modelRef\n";
	      }
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

  my $coordinations = {};
  my $shells = $self->{shells};
  foreach my $shell (@$shells)
    {
    my @models;
    foreach my $major ("Tetrahedral", "TrigonalBipyramidal", "Octahedral", "PentagonalBipyramidal")
      {
      my $cgObj = $major->new(shellObj => $shell);
      $cgObj->bestTestStatistic("dev");
      push @models, $cgObj;
      }

    my $bestModel = (sort {$a->{bestCombo}->{deviation} <=> $b->{bestCombo}->{deviation}} (grep {defined $_->{bestCombo} && $_->{bestCombo}->{deviation} != 0;} (@models)))[0];
    next if (! $bestModel);
    my $modelRef = ref $bestModel;
    push @{$$coordinations{$modelRef}}, $bestModel;
    }

  $self->{coordinations} = $coordinations;
  }


## print out the counts of each CGs on screen
sub printStats
  {
  my $self = shift @_;
  my $i = (@_)? shift @_ : 0;

  my $now_string = localtime;
  print "$now_string\n";
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

  my $coordinationAngles = {};
  foreach my $coord (keys %$coordSets)
    {
    my $coordSet = $$coordSets{$coord};
    next if (! $coordSet);
    map { $_->calcExpectedAngleStats($coordinationAngles) ;} (@$coordSet);
    }

  my $totalDev = 0;
  my $totaln = 0;
  foreach my $coordination (keys %$coordinationAngles)
    {
    foreach my $angle (keys %{$$coordinationAngles{$coordination}})
      {
      my $stats = RawStatistics->new("variables" => $$coordinationAngles{$coordination}{$angle} ) ;
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
      $$angleStats{$coordination}{$angle}{"variance"} = $stats->calcVariance;
      $$angleStats{$coordination}{$angle}{"standardDeviation"} = $stats->calcStd();
      #$$angleStats{$coordination}{$angle}{"min"} = $stats->min();
      #$$angleStats{$coordination}{$angle}{"max"} = $stats->max();
      }
    }

  my ($pooledVar, $pooledCount);
  map {my $cg = $_; map {$pooledCount += $$angleStats{$cg}{$_}{"count"} - 1; $pooledVar += $$angleStats{$cg}{$_}{"variance"} * ($$angleStats{$cg}{$_}{"count"} - 1) } (keys %{$$coordinationAngles{$cg}})} (keys %$coordinationAngles);
  $pooledVar = $pooledVar / $pooledCount;

  my $variance = $totalDev/$totaln ;
  $$angleStats{"variance"} = $variance;
  $$angleStats{"standardDeviation"} = $variance ** 0.5;
  $$angleStats{"pooledVar"} = $pooledVar;
  $$angleStats{"pooledStd"} = $pooledVar ** 0.5;

  $self->{rawAngles} = $coordinationAngles;
  return $angleStats;
  }

## Given coordination classification, calculated distance statistics       
sub calcDistStats
  {
  my $self = shift @_;

  my $coordSets = $self->{coordinations};
  my $elementDists = {};
  my $distStats = {};
  my $elementRes = {};

  foreach my $coord (keys %$coordSets)
    {
    my $coordSet = $$coordSets{$coord};

    foreach my $model (@$coordSet)
      {
      foreach my $ligand ( @{$model->{bestCombo}->{ligands}} )
        {
        my $element = $ligand->{element};
        my $resolution = ($ligand->{resolution} == -1)? 2.5 : $ligand->{resolution};

        push (@{$$elementDists{$element}}, $model->{shellObj}->{center}->distance($ligand));
        push (@{$$elementRes{$element}}, $resolution);
        }
      }
    }

  foreach my $element (keys %$elementDists)
    {
    map {push (@{$$elementDists{"average"}}, $_) ;} (@{$$elementDists{$element}});
    map {push (@{$$elementRes{"average"}}, $_) ;} (@{$$elementRes{$element}});

    if (@{$$elementDists{$element}} > 30)
      {
      my $stats = RawStatistics->new("variables" => $$elementDists{$element} ) ;
      my $statsRes = RawStatistics->new("variables" => $$elementRes{$element} ) ;

      $$distStats{$element} = {"mean" => $stats->calcMean(), "variance" => $stats->calcVariance(), "count" => $stats->count(), "max" => $stats->max(), "min" => $stats->min(), "standardDeviation" => $stats->calcStd(), "resolutionAvg" => $statsRes->calcMean()};
      }
    }

  my ($pooledVar, $pooledCount);
  map {$pooledCount += $$distStats{$_}{"count"} - 1; $pooledVar += $$distStats{$_}{"variance"} * ($$distStats{$_}{"count"} - 1) } (keys %$elementDists);
  $pooledVar = $pooledVar / $pooledCount;

  my $stats = RawStatistics->new("variables" => $$elementDists{"average"} ) ;
  my $statsRes = RawStatistics->new("variables" => $$elementRes{"average"} ) ;

  $$distStats{"average"} = {"mean" => $stats->calcMean(), "variance" => $stats->calcVariance(), "count" => $stats->count(), "max" => $stats->max(), "min" => $stats->min(), "standardDeviation" => $stats->calcStd(), "resolutionAvg" => $statsRes->calcMean(), "pooledVar" => $pooledVar };

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












 
