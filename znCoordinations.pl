#!/usr/bin/perl 

my $help = <<HELP;
  usage: ./znCoordinations.pl PDBpathsFile -i/d [(criteria threshold) statisticsFilePrefix] [-s sequences header sequenceResultsFile] [-r rSpreadsheet (leaveOut) statisticsFile]
		[-json jsonFile] [-dumper dumperFile] [-decision] [-angleList angleListMid] [-ec pathToFlat pathToPDB ecFile] [-ligand] [-angleBD angleBreakDownFile]
 
  Parameters:

    required:
  	PDBpathsFile 
	-i/d [(criteria threshold) statisticsFilePrefix], 
		criteria: probability(p)/compressed(c)/nonModel(n); 
		threshold: threshold for the criteria; 
		statisticsFilePrefix: the prefix of the statistics file

    optional:
	[-s sequences header sequenceResultsFile], sequence: SEQRES(s)/ATOM(a)/numbering(n); header: binding(b)/shell(s)/non-aa to second-shell aa(ss)/non-aa to closet aa(c)
	[-r rSpreadsheet (leaveOut) statisticsFile], rSpreadsheet: the file for R input; (leaveOut): leaveOut/leave/l, whether to leave out the smallest angle in the chi-squared prob cals; statisticsFile: the statistics file for chi-squared prob calculation
	[-json jsonFile], this is an option only for after nonModel (-i/d n threshold prefix)
	[-dumper dumperFile], this is an option only for after nonModel (-i/d n threshold prefix)
	[-decision], the output will be print on screen, and it is an option only for after IA (-i)
	[-angleList angleListMid], angleListMid: the name root for
	[-ec pathToFlat pathToPDB ecFile], pathToFlat: path to uniprot flat file; pathToPDB, path to PDB parent diredtory; ecFile: the output file
	[-ligand], the output will be print on screen
	[-angleBD angleBreakDownFile]
HELP

#===============================================================================
#
#         FILE:  "znCoordinations.pl"
#
#        USAGE:  ./znCoordinations.pl
#
#  DESCRIPTION: 
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Sen Yao 
#      COMPANY:  
#      VERSION:  0.1 
#      CREATED:  1/15/2013 02:07:30 PM
#     REVISION:  ---
#===============================================================================

use strict;
use BasicTable2 qw(:ALL);
use ZnCGanalysis;

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

use JSON;
use Data::Dumper::Concise;

########## read in required arguments ##########
my $pathsFile = shift @ARGV;
unless ($ARGV[0] eq "-i" || $ARGV[0] eq "-d")
  { print STDERR $help; exit;}
my $flow = shift @ARGV;

## This is not really -i/-d specific
my ($iaControl, $threshold); 
if ($ARGV[0] =~ /^(p|porb|probability)$/)
  {
  shift @ARGV;
  $iaControl = "p";
  if ($ARGV[0] =~ /^[+-]?\d+(\.\d+)?$/)
    { $threshold = shift @ARGV; }
  else 
    { $threshold = 0.001;}
  }
elsif ($ARGV[0] =~ /^(c|comp|compressed)$/)
  {
  shift @ARGV;
  $iaControl = "c";
  if ($ARGV[0] =~ /^[+-]?\d+(\.\d+)?$/)
    { $threshold = shift @ARGV; }
  else
    { $threshold = 68;}
  }
## nonModel is to remove the zinc shell for all later iterations if if didn't pass the probability threshold.
elsif ($ARGV[0] =~ /^(n|non|nonModel)$/) 
  { 
  shift @ARGV;
  $iaControl = "n";
  if ($ARGV[0] =~ /^[+-]?\d+(\.\d+)?$/)
    { $threshold = shift @ARGV; }
  else
    { $threshold = 0.001;}
  }

## Assign names for statistics file for each round
my $statOutFileName = ($ARGV[0] !~ /\-/) ? shift @ARGV : "statistics";

########## read in optional arguments ##########
my ($seqOpt, $seq, $header, $seqFile, $rInputOpt, $rInputFile, $statsFile, $leaveOut, $jsonOpt, $jsonFile, $dumperOpt, $dumperFile, $decisionOpt, $angleListOpt, $angleListMid);
my ($ECopt, $pathToFlat, $pathToPDB, $ecFile, $ligandOpt, $angleBreakDownOpt, $angleBreakDownFile);
while (@ARGV) 
  {
  if ($ARGV[0] !~ /\-/)
    { print STDERR $help; exit;}
  my $option = shift @ARGV;

  if ($option eq "-s")
    {
    $seqOpt = 1;
    $seq = shift @ARGV;
    $header = shift @ARGV;
    $seqFile = shift @ARGV;
    
    unless ($seq eq "s" || $seq eq "a" || $seq eq "n")
      { print STDERR $help; exit;}
    unless ($header eq "b" || $header eq "s" || $header eq "ss" || $header eq "c")
      { print STDERR $help; exit;} 
    }
  elsif ($option eq "-r")
    {
    $rInputOpt = 1;
    $rInputFile = shift @ARGV;
    $leaveOut = shift @ARGV if ($ARGV[0] =~ /^(l|leave|leaveOut)$/);
    $statsFile = shift @ARGV if ($ARGV[0] !~ /\-/);
    }
  elsif ($option eq "-json")
    {
    $jsonOpt = 1;
    $jsonFile = shift @ARGV;
    }
  elsif ($option eq "-dumper")
    {
    $dumperOpt = 1;
    $dumperFile = shift @ARGV;  
    }
  elsif ($option eq "-decision")
    {
    $decisionOpt = 1;
    }
  elsif ($option eq "-angleList")
    {
    $angleListOpt = 1;
    $angleListMid = shift @ARGV;
    }
  elsif ($option eq "-ec")
    {
    $ECopt = 1;
    $pathToFlat = shift @ARGV;  # Path to the file of all unpid and ec number generated from uniprot flat file
    $pathToPDB = shift @ARGV;  # with or without '/' at the end
    $ecFile = shift @ARGV;
    }
  elsif ($option eq "-ligand")
    {
    $ligandOpt = 1;
    }
  elsif ($option eq "-angleBD")
    {
    $angleBreakDownOpt = 1;
    $angleBreakDownFile = shift @ARGV;
    }
  }


######################################## 
############# main process #############
########################################

my $coord = ZnCGanalysis->new("pathsFile" => $pathsFile, "element" => "ZN");
print "Zn: ", $coord->{numCenter}, "\n";
print "Cluster: ", $coord->{numCluster}, "\n";
$coord->bootstrapCoordination($statOutFileName);

## required argument defines the workflow
if ($flow eq "-i") 
  { $coord->IAcoordination($statOutFileName, $iaControl, $threshold); }
else 
  { $coord->bindShellViaDist($statOutFileName); }

## optional args set what to print out
if ($seqOpt)
  {
  $coord->printSequences($seqFile, $seq, $header, "four"); ## four ligands
  }

if ($rInputOpt)
  {
  open (FID, ">", $rInputFile) or die $!;
  my $stats = &readTableFile($statsFile) if ($statsFile);

  print FID "Zn_ID\tAngle_1\tAngle_2\tAngle_3\tAngle_4\tAngle_5\tAngle_6\tLigand_1\tLigand_2\tLigand_3\tLigand_4\t";
  print FID "Bond_Length_1\tBond_length_2\tBond_length_3\tBond_length_4\tBi_status_1\tBi_status_2\tBi_status_3\tBi_status_4\tBi_status_5\tBi_status_6\t";
  print FID "Method\tDate\tResolution\tB_factor_1\tB_factor_2\tB_factor_3\tB_factor_4\tOverall_bi_status\tLig_chain_ID_combo\tLig_residue_combo\tLig_atom_combo\tAmineN\t";
  print FID "Chi_prob_Tet\tChi_prob_Bva\tChi_prob_Bvp\tChi_prob_Spv\tChi_prob_Spl\t" if ($statsFile);
  print FID "\n";

  foreach my $znObj (@{$coord->{coordinations}{"four"}})
    {
    print FID $znObj->{shellObj}->znID(), "\t";
    map { print FID "$_  \t";} ($znObj->orderedAngles());
    map { print FID "$_  \t";} ($znObj->ligandAtomElement()); ## also have bond lengths for each ligand and bidentations for each angle
    
    print FID $znObj->{shellObj}->{center}->{method}, "\t";
    print FID $znObj->{shellObj}->{center}->{date}, "\t";
    print FID $znObj->{shellObj}->{center}->{resolution}, "\t";
    map { print FID $_->{bFactor}, "\t";} (@{$znObj->{bestCombo}->{ligands}});
    print FID $znObj->bidentate(), "\t";
    map { print FID "$_  \t";} ($znObj->ligandCombos());

    map { print FID $_, "\t";} (&coordProbs($znObj->{shellObj}, $stats, $leaveOut)) if ($statsFile);
    #map { print FID $_->{chiAngle}, "\t";} (@{$znObj->{bestCombo}->{ligands}});
    
    print FID "\n";
    }
  close FID;
  }

sub coordProbs
  {
  my $shell = shift @_;
  my $stats = shift @_;
  my $leaveOut = shift @_;

  my $th = Tetrahedral->new(shellObj => $shell);
  my $bva = TrigonalBipyramidalVA->new(shellObj => $shell);
  my $bvp = TrigonalBipyramidalVP->new(shellObj => $shell);
  my $spv = SquarePyramidalV->new(shellObj => $shell);
  my $spl = SquarePlanar->new(shellObj => $shell);

  $th->bestTestStatistic("ownStats", 0, 0, $leaveOut, $stats);
  $bva->bestTestStatistic("ownStats", 0, 0, $leaveOut, $stats);
  $bvp->bestTestStatistic("ownStats", 0, 0, $leaveOut, $stats);
  $spv->bestTestStatistic("ownStats", 0, 0, $leaveOut, $stats);
  $spl->bestTestStatistic("ownStats", 0, 0, $leaveOut, $stats);

  my @probs;
  push @probs, $th->{bestCombo}->{probability};
  push @probs, $bva->{bestCombo}->{probability};
  push @probs, $bvp->{bestCombo}->{probability};
  push @probs, $spv->{bestCombo}->{probability};
  push @probs, $spl->{bestCombo}->{probability};

  return @probs;
  }


## Some other supplemental options
if ($jsonOpt)
  {
  open (JOUT, '>', $jsonFile);
  my $jsonObj = JSON->new->allow_blessed->convert_blessed->encode( $coord->{nonModels} );
  print JOUT $jsonObj;
  }

if ($dumperOpt)
  {
  open (DOUT, '>', $dumperFile );
  print DOUT Dumper($coord->{nonModels});
  }

if ($decisionOpt)
  {
  foreach my $key ( sort keys %{$coord->{decisions}} )
    { print "$key, ", $coord->{decisions}{$key}, "\n"; }
  }

if ($angleListOpt)
  {
  foreach my $coordination (keys %{$coord->{coordinations}} )
    {
    next if (! $coord->{coordinations}{$coordination});

    my $fileName = "$coordination.$angleListMid.txt";
    open (FID, ">", $fileName) or die $!;
    foreach my $znObj (@{$coord->{coordinations}{$coordination}})
      {
      map { print FID  "$_  \t";} ($znObj->angleList());
      my $prob = $znObj->{bestCombo}->{probability};
      print FID "$prob\n";
      print FID "\n";
      }
    close FID;
    }
  }

if ($ECopt)
  {
  die if (! open (UNPFLAT, $pathToFlat));
  $pathToPDB = $1 if ($pathToPDB =~ /^(.+)\/$/);

  my $ecUnp = {};
  foreach my $coordination (keys %{$coord->{coordinations}})
    {
    foreach my $znObj (@{$coord->{coordinations}{$coordination}})
      {
      my $model = ref $znObj;
      my $znAtom = $znObj->{shellObj}->{center};
      my $pdbid = $znAtom->{PDBid};
      my $chainid = $znAtom->{chainID};
      my $serial = $znAtom->{residueNumber};

      my $id = "$pdbid.$chainid.$serial";
      my $inputId = lc($pdbid).uc($chainid);
  
      $$ecUnp{$id}{"uniprotID"} = [&getUnpidFromPdbid($inputId)];
      $$ecUnp{$id}{"coordination"} = $model;
      }  
    }

  while(my $input = <UNPFLAT>) #foreach uniprot id, check if it's zn-associated
    {
    chomp $input;
    next if ($input =~ /^\s*$/);
    next if ($input !~ /,/); #if the uniprot doesn't have ec # at all

    my %unpid = map {($_, 1);} (split (/; /, $input));
    foreach my $id (keys %$ecUnp)
      {
      if ( grep {$unpid{$_};} (@{$$ecUnp{$id}{"uniprotID"}}) ) # accession number matches => it's zn-associated
        {
        my @numbers = split (/, /, $input);
        foreach my $number (@numbers)
          {
          my $oneEC = {};

          if ($number =~ /(.)\.(.+)\.(.+)\.(.+)$/)
            {
            $$oneEC{"id"} = $1.".".$2.".".$3.".".$4;
            next if grep {$$_{"id"} eq $$oneEC{"id"};} (@{$$ecUnp{$id}{"ecNumbers"}});

            $$oneEC{"level1"} = $1;
            $$oneEC{"level2"} = $2;
            $$oneEC{"level3"} = $3;
            $$oneEC{"level4"} = $4;
            $$oneEC{"from"} = "unpflat";

            push (@{$$ecUnp{$id}{"ecNumbers"}}, $oneEC);
	    }
          }
        }
      }
    }
  
  &writeTableFile($ecFile, $ecUnp);
  }

sub openPdbfile
  {
  my $pdbid = shift @_;
  my ($filename, $foldername, $chain) ;

  if ( $pdbid =~ /^(\w(\w\w)\w)(\w)$/ )
    {
    $filename = $1;
    $foldername = $2;
    $chain = $3;
    }
  else
    {print "error in pdb id format: $pdbid\n";}
  
  my $path = "$pathToPDB/$foldername/pdb$filename\.ent\.gz";
  die if (! open (my $pdbfile, "zcat $path|"));
  return ($pdbfile, $filename, $chain);
  }


sub getUnpidFromPdbid
  {
  my $pdbid = shift @_;
  my ($pdbfile, $pdbid, $chain) = &openPdbfile($pdbid);

  my @uniprots;
  my %unpid;
  while (my $line = <$pdbfile>)
    {
    chomp $line;
    if ($line =~ /^DBREF.+ $chain .+ UNP/)
      {
      my $unp = substr ($line, 33, 6);
      $unpid{$unp} = 1;
      }
    }
  foreach my $id (keys %unpid)
    {push @uniprots, $id;}
  return @uniprots;
  }

if ($ligandOpt)
  {
  my %ligandCombo;
  foreach my $coordination (values %{$coord->{coordinations}})
    {
    foreach my $znObj (@$coordination)
      {
      my @residues = map {$_->{residueName};} (@{$znObj->{bestCombo}->{ligands}});
      $ligandCombo{join ".", sort {lc($a) cmp lc($b);} (@residues)}++;
      }
    }
  map {print "$_= ", $ligandCombo{$_}, "\n";} (sort {$ligandCombo{$b} <=> $ligandCombo{$a}} (keys %ligandCombo));
  }

if ($angleBreakDownOpt)
  {
  foreach my $coordination (keys %{$coord->{rawAngles}})
    {
    foreach my $angel (keys %{$coord->{rawAngles}{$coordination}} )
      {
      my $histogram=[];

      push @$histogram, scalar (grep {$_ > 0  && $_ < 15;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 15 && $_ < 25;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 25 && $_ < 35;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 35 && $_ < 45;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 45 && $_ < 55;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 55 && $_ < 65;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 65 && $_ < 75;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 75 && $_ < 85;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 85 && $_ < 95;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 95 && $_ < 105;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 105 && $_ < 115;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 115 && $_ < 125;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 125 && $_ < 135;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 135 && $_ < 145;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 145 && $_ < 155;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 155 && $_ < 165;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 165 && $_ < 175;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 175 && $_ < 180;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 180 ;} (@{$coord->{rawAngles}{$coordination}{$angel}}) ) ;
 
      $coord->{rawAngles}{$coordination}{$angel} = $histogram;
      }
    }

  &writeTableFile($angleBreakDownFile, $coord->{rawAngles});
  }




