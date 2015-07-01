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

  my ($three, $four, $five, $six, $more);
  foreach my $shell (@{$self->{shells}})
    {
      my $tpl = TrigonalPlanar->new(shellObj => $shell);
      my $th = Tetrahedral->new(shellObj => $shell);
      my $tb = TrigonalBipyramidal->new(shellObj => $shell);
      my $oh = Octahedral->new(shellObj => $shell);

      $tpl->bestDistChi($stats);
      $th->bestDistChi($stats);
      $tb->bestDistChi($stats);
      $oh->bestDistChi($stats);

    if (! defined $tpl->{bestCombo}) ## less than 3 ligands
      {
      next;
      }
    elsif (! defined $th->{bestCombo}) ## 3 ligands
      {
      push @$three, $tpl;
      }
    elsif (! defined $tb->{bestCombo}) ## 4 ligands
      {
      push @$four, $th;
      }
    elsif (! defined $oh->{bestCombo}) ## 5 ligands
      {
      my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($th, $tb)))[0];
      if (ref $bestModel eq "Tetrahedral")
        { push @$four, $th;}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { push @$five, $tb;}
      }
    else
      {
      my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($th, $tb, $oh)))[0];
      if (ref $bestModel eq "Tetrahedral")
        { push @$four, $th;}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { push @$five, $tb;}
      elsif (ref $bestModel eq "Octahedral")
        { push @$six, $oh;}
      }
 
    ## print chi probabilities
    #print $shell->znID(), ": ";
    #print "th, ", $th->{bestCombo}->{probability}, "; ";
    #print "tb, ", $tb->{bestCombo}->{probability}, "; ";
    #print "oh, ", $oh->{bestCombo}->{probability}, "; ";
    #print "\n";
    }

  &writeTableFile("$statOutFileName.dist.txt", $stats);

  my $bestDist = {};
  %$bestDist = (  "three" => $three,
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
    foreach my $lig (@{$model->{bestCombo}->{ligands}})
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

## Classify coordination using chi statistics
sub calcChiCoordination
  {
  my $self = shift @_;
  my $control = shift @_;
  my $threshold = shift @_;
  my $stats = (@_)? (shift @_) : ($self->{stats});

  my ($tetrahedrals, $trigonalbipyramidals, $octahedrals, $trigonalplanars, $tetrahedralVs, $trigonalbipyramidalVAs, $trigonalbipyramidalVPs, $squareplanars, $squarepyramidalVs, $squarepyramidals, $unusables);
  my $decisions = {};
  foreach my $shell (@{$self->{shells}})
    {
    # major coordinations
    my $th = Tetrahedral->new(shellObj => $shell);
    my $tb = TrigonalBipyramidal->new(shellObj => $shell);
    my $oh = Octahedral->new(shellObj => $shell);

    # 3 ligands
    my $tpl = TrigonalPlanar->new(shellObj => $shell);
    my $thv = TetrahedralV->new(shellObj => $shell);

    # 4 ligands, trignal bipyramidal related
    my $bva = TrigonalBipyramidalVA->new(shellObj => $shell);
    my $bvp = TrigonalBipyramidalVP->new(shellObj => $shell);

    # 4 ligands, square pyramidal related
    my $spv = SquarePyramidalV->new(shellObj => $shell);
    my $spl = SquarePlanar->new(shellObj => $shell);

    # 5 ligands
    my $spy = SquarePyramidal->new(shellObj => $shell);

    $th->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $tb->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $oh->bestTestStatistic("chi", $control, $threshold, 0, $stats);

    $tpl->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $thv->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $bva->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $bvp->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spv->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spl->bestTestStatistic("chi", $control, $threshold, 0, $stats);
    $spy->bestTestStatistic("chi", $control, $threshold, 0, $stats);


    if (! defined $tpl->{bestCombo}) ## less than 3 ligands
      {
      $$decisions{"012"}++;
      next;
      }
    elsif (! defined $th->{bestCombo}) ## 3 ligands
      {
      my $bestModel = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($tpl, $thv)))[0];

      if ($bestModel->{bestCombo}->{probability} > 0.01)
	{
	if (ref $bestModel eq "TetrahedralV" && $bestModel)
          {
          push @$tetrahedralVs, $thv;
          $$decisions{"3Thv"} += 1;
          }
        elsif (ref $bestModel eq "TrigonalPlanar")
          {
          push (@$trigonalplanars, $tpl) ;
          $$decisions{"3Tpl"} += 1;
  	  }
 	}
      else
	{$$decisions{"3None"} += 1;}
      next;
      }
    elsif (! defined $tb->{bestCombo}) ## 4 ligands
      {
      my @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($th, $bva, $bvp, $spv, $spl)));
      my $bestModel = $models[0] ;
      my $sceondModel = $models[1] ;
      my $thirdModel = $models[2] ;


      ## set the probability threshold, remove low prob ones from statistics calculation. 
      if ($control eq "p" && $bestModel->{bestCombo}->{probability} < $threshold) 
	{
	push @$unusables, $bestModel;
        $$decisions{"4Unusable"} += 1;
	next;
	}

      if (ref $bestModel eq "Tetrahedral")
        { 
	push @$tetrahedrals, $th; 
	$$decisions{"4Th"} += 1;
	}
      elsif (ref $bestModel eq "TrigonalBipyramidalVA")
        {
	if (ref $sceondModel eq "Tetrahedral" && ($bestModel->{bestCombo}->{probability} < (2 * $sceondModel->{bestCombo}->{probability})) )
	  {
          push @$tetrahedrals, $th;
          $$decisions{"4Th"} += 1;
          }
	else
 	  {
 	  push @$trigonalbipyramidalVAs, $bva; 
	  $$decisions{"4Bva"} += 1;
	  }
	}
      elsif (ref $bestModel eq "TrigonalBipyramidalVP")
        {
        push @$trigonalbipyramidalVPs, $bvp;
        $$decisions{"4Bvp"} += 1;
        }
      elsif (ref $bestModel eq "SquarePyramidalV.pm")
        {
        push @$squarepyramidalVs, $spv;
        $$decisions{"4Spv"} += 1;
        }
      elsif (ref $bestModel eq "SquarePlanar")
        { 
	push @$squareplanars, $spl; 
	$$decisions{"4Spl"} += 1;
	}

      next;
      }
    elsif (! defined $oh->{bestCombo}) ## 5 ligands
      {
      my @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($th, $tb, $bva, $bvp, $spv, $spl, $spy)));
      my $bestModel = $models[0] ;
      my $sceondModel = $models[1] ;
      my $thirdModel = $models[2] ;

      ## Save as 4 ligands above
      if ($control eq "p" && $bestModel->{bestCombo}->{probability} < $threshold)
        {
        push @$unusables, $bestModel;
        $$decisions{"5Unusable"} += 1;
	next;
        }

      if (ref $bestModel eq "Tetrahedral")
        { 
	push @$tetrahedrals, $th; 
	$$decisions{"5Th"} += 1;
	}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { 
	push @$trigonalbipyramidals, $tb; 
	$$decisions{"5Tb"} += 1;
	}

      # trigonal bipyramidal related
      elsif (ref $bestModel eq "TrigonalBipyramidalVA")
        {
        if (ref $sceondModel eq "TrigonalBipyramidal") 
	  { 
  	  push @$trigonalbipyramidals, $tb; 
	  $$decisions{"5BvaTb"} += 1;
	  }
        elsif (ref $sceondModel eq "TrigonalBipyramidalVP" && ref $thirdModel eq "TrigonalBipyramidal" )
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"5BvaBvpTb"} += 1;
          }
	elsif (ref $sceondModel eq "Tetrahedral" && ($bestModel->{bestCombo}->{probability} < (2 * $sceondModel->{bestCombo}->{probability})) )
          {
          push @$tetrahedrals, $th;
          $$decisions{"5Th"} += 1;
          }
        else
	  { 
	  push @$trigonalbipyramidalVAs, $bva; 
	  $$decisions{"5Bva"} += 1;
	  }
	}
      elsif (ref $bestModel eq "TrigonalBipyramidalVP")
        {
        if (ref $sceondModel eq "TrigonalBipyramidal")
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"5BvpTb"} += 1;
          }
        elsif (ref $sceondModel eq "TrigonalBipyramidalVA" && ref $thirdModel eq "TrigonalBipyramidal" )
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"5BvpBvaTb"} += 1;
          }
        else
          {
          push @$trigonalbipyramidalVPs, $bvp;
          $$decisions{"5Bvp"} += 1;
          }
        }

      # square pyramidal related
      elsif (ref $bestModel eq "SquarePyramidal")
        {
        push @$squarepyramidals, $spy;
        $$decisions{"5Spy"} += 1;
        }
      elsif (ref $bestModel eq "SquarePyramidalV")
        {
        if (ref $sceondModel eq "SquarePyramidal")
          {
          push @$squarepyramidals, $spy;
          $$decisions{"5SpvSpy"} += 1;
          }
        elsif (ref $sceondModel eq "SquarePlanar" && ref $thirdModel eq "SquarePyramidal" )
          {
          push @$squarepyramidals, $spy;
          $$decisions{"5SpvSplSpy"} += 1;
	  }
        else
          {
          push @$squarepyramidalVs, $spv;
          $$decisions{"5Spv"} += 1;
          }
        }
      elsif (ref $bestModel eq "SquarePlanar")
	{
        if (ref $sceondModel eq "SquarePyramidal")
          { 
	  push @$squarepyramidals, $spy;
	  $$decisions{"5SplSpy"} += 1;
	  }
        elsif (ref $sceondModel eq "SquarePyramidalV" && ref $thirdModel eq "SquarePyramidal" )
          {
          push @$squarepyramidals, $spy;
          $$decisions{"5SplSpvSpy"} += 1;
          }
        else
          { 
	  push @$squareplanars, $spl; 
	  $$decisions{"5Spl"} += 1;
	  }
	}

      next;
      }

    else ## 6 ligands or more
      { 
      my @models = (sort {$b->{bestCombo}->{probability} <=> $a->{bestCombo}->{probability}} (grep {$_->{bestCombo}->{probability} != 0;} ($th, $tb, $oh, $bva, $bvp, $spv, $spl, $spy)));
      my $bestModel = $models[0] ;
      my $sceondModel = $models[1] ;
      my $thirdModel = $models[2];
      my $fourthModel = $models[3];

      ## Same as above
      if ($control eq "p" && $bestModel->{bestCombo}->{probability} < $threshold)
        {
        push @$unusables, $bestModel;
        $$decisions{"6Unusable"} += 1;
        next;
	}

      if (ref $bestModel eq "Tetrahedral")
        { 
	push @$tetrahedrals, $th; 
	$$decisions{"6Th"} += 1;
	}
      elsif (ref $bestModel eq "TrigonalBipyramidal")
        { 
	push @$trigonalbipyramidals, $tb; 
	$$decisions{"6Tb"} += 1;
	}
      elsif (ref $bestModel eq "Octahedral")
        { 
	push @$octahedrals, $oh; 
	$$decisions{"6Oh"} += 1;
	}

      # trigonal bipyramidal related
      elsif (ref $bestModel eq "TrigonalBipyramidalVA")
        {
        if (ref $sceondModel eq "TrigonalBipyramidal")
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"6BvaTb"} += 1;
          }
        elsif (ref $sceondModel eq "TrigonalBipyramidalVP" && ref $thirdModel eq "TrigonalBipyramidal" )
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"6BvaBvpTb"} += 1;
          }
        elsif (ref $sceondModel eq "Tetrahedral" && ($bestModel->{bestCombo}->{probability} < (2 * $sceondModel->{bestCombo}->{probability})) )
          {
          push @$tetrahedrals, $th;
          $$decisions{"6Th"} += 1;
          }
        else
          {
          push @$trigonalbipyramidalVAs, $bva;
          $$decisions{"6Bva"} += 1;
          }
        }
      elsif (ref $bestModel eq "TrigonalBipyramidalVP")
        {
        if (ref $sceondModel eq "TrigonalBipyramidal")
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"6BvpTb"} += 1;
          }
        elsif (ref $sceondModel eq "TrigonalBipyramidalVA" && ref $thirdModel eq "TrigonalBipyramidal" )
          {
          push @$trigonalbipyramidals, $tb;
          $$decisions{"6BvpBvaTb"} += 1;
          }
        else
          {
          push @$trigonalbipyramidalVPs, $bvp;
          $$decisions{"6Bvp"} += 1;
          }
        }

      # octahedral related
      elsif (ref $bestModel eq "SquarePyramidal")
        {
        if (ref $sceondModel eq "Octahedral")
          {
          push @$octahedrals, $oh;
          $$decisions{"6SpyOh"} += 1;
          }
        else
          {
          push @$squarepyramidals, $spy;
          $$decisions{"6Spy"} += 1;
          }
        }
      elsif (ref $bestModel eq "SquarePyramidalV")
        {
        if (ref $sceondModel eq "Octahedral")
          {
          push @$octahedrals, $oh;
          $$decisions{"6SpvOh"} += 1;
          }
        elsif (ref $sceondModel eq "SquarePyramidal")
          {
          if (ref $thirdModel eq "Octahedral")
            {
            push @$octahedrals, $oh;
            $$decisions{"6SpvSpyOh"} += 1;
            }
          else
            {
            push @$squarepyramidals, $spy;
            $$decisions{"6SpvSpy"} += 1;
            }
          }
        elsif (ref $sceondModel eq "SquarePlanar" && ref $thirdModel eq "SquarePyramidal" )
          {
          if (ref $fourthModel eq "Octahedral")
            {
            push @$octahedrals, $oh;
            $$decisions{"6SpvSplSpyOh"} += 1;
            }
          else
            {
            push @$squarepyramidals, $spy;
            $$decisions{"6SpvSplSpy"} += 1;
            }
	  }
        else
          {
          push @$squarepyramidalVs, $spv;
          $$decisions{"6Spv"} += 1;
          }
        }
      elsif (ref $bestModel eq "SquarePlanar")
        {
        if (ref $sceondModel eq "Octahedral")
          {
          push @$octahedrals, $oh;
          $$decisions{"6SplOh"} += 1;
          }
        elsif (ref $sceondModel eq "SquarePyramidal")
          {
          if (ref $thirdModel eq "Octahedral")
            {
            push @$octahedrals, $oh;
            $$decisions{"6SplSpyOh"} += 1;
            }
          else
            {
            push @$squarepyramidals, $spy;
            $$decisions{"6SplSpy"} += 1;
            }
          }
        elsif (ref $sceondModel eq "SquarePyramidalV" && ref $thirdModel eq "SquarePyramidal" )
          {
          if (ref $fourthModel eq "Octahedral")
            {
            push @$octahedrals, $oh;
            $$decisions{"6SplSpvSpyOh"} += 1;
            }
          else
            {
            push @$squarepyramidals, $spy;
            $$decisions{"6SplSpvSpy"} += 1;
            }
          }
        else
          {
          push @$squareplanars, $spl;
          $$decisions{"6Spl"} += 1;
          }
        }

      next;
      }
    }

  my $coordinations = {};
  %$coordinations = (   "tetrahedrals" => $tetrahedrals,
			"trigonalBipyramidals" => $trigonalbipyramidals,
			"octahedrals" => $octahedrals, 
			"trigonalPlanars" => $trigonalplanars,
                        "tetrahedralsVacancy" => $tetrahedralVs,
                        "trigonalBipyramidalsVacancyAxial" => $trigonalbipyramidalVAs,
                        "trigonalBipyramidalsVacancyPlanar" => $trigonalbipyramidalVPs,
                        "squarePyramidalsVacancy" => $squarepyramidalVs,
			"squarePlanars" => $squareplanars,
			"squarePyramidals" => $squarepyramidals );

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
    my $th = Tetrahedral->new(shellObj => $shell);
    my $tb = TrigonalBipyramidal->new(shellObj => $shell);
    my $oh = Octahedral->new(shellObj => $shell);

    $th->bestTestStatistic("dev");
    $tb->bestTestStatistic("dev");
    $oh->bestTestStatistic("dev");

    my $bestModel = (sort {$a->{bestCombo}->{deviation} <=> $b->{bestCombo}->{deviation}} (grep {$_->{bestCombo}->{deviation} != 0;} ($th, $tb, $oh)))[0];
    if (ref $bestModel eq "Tetrahedral")
      { push @$tetrahedrals, $th; }
    elsif (ref $bestModel eq "TrigonalBipyramidal")
      { push @$trigonalbipyramidals, $tb; }
    elsif (ref $bestModel eq "Octahedral")
      { push @$octahedrals, $oh; }

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

#  print $self->{coordinations}->{tetrahedrals}, "flag\n";
  print "\nTetrahedral = ", scalar @{$self->{coordinations}->{tetrahedrals}} if $self->{coordinations}->{tetrahedrals}; 
  print "\nTrigonal bipyramidal = ", scalar @{$self->{coordinations}->{trigonalBipyramidals}} if  $self->{coordinations}->{trigonalBipyramidals};
  print "\nOctahedral = ", scalar @{$self->{coordinations}->{octahedrals}} if $self->{coordinations}->{octahedrals};
  print "\nTrigonal planar = ", scalar @{$self->{coordinations}->{trigonalPlanars}} if $self->{coordinations}->{trigonalPlanars};
  print "\nTetrahedral with vacancy = ", scalar @{$self->{coordinations}->{tetrahedralsVacancy}} if $self->{coordinations}->{tetrahedralsVacancy};
  print "\nTrigonal bipyramidal with vacancy axial = ", scalar @{$self->{coordinations}->{trigonalBipyramidalsVacancyAxial}} if  $self->{coordinations}->{trigonalBipyramidalsVacancyAxial};
  print "\nTrigonal bipyramidal with vanancy planar = ", scalar @{$self->{coordinations}->{trigonalBipyramidalsVacancyPlanar}} if  $self->{coordinations}->{trigonalBipyramidalsVacancyPlanar};
  print "\nSquare pyramidal with vancy = ", scalar @{$self->{coordinations}->{squarePyramidalsVacancy}} if $self->{coordinations}->{squarePyramidalsVacancy};
  print "\nSquare planars = ", scalar @{$self->{coordinations}->{squarePlanars}} if  $self->{coordinations}->{squarePlanars};
  print "\nSquare pyramidal = ", scalar @{$self->{coordinations}->{squarePyramidals}} if $self->{coordinations}->{squarePyramidals};
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












 
