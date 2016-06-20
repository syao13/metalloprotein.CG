#!/usr/bin/perl 

my $help = <<HELP;
  usage: ./metalCoordinations.pl PDBpathsFile metal [shellCutoff shellElmt] -i/d/dd/bs [(criteria threshold) statisticsFilePrefix] 
		[-s sequences header sequenceResultsFile] [-r rSpreadsheet (leaveOut) statisticsFile]
		[-json jsonFile] [-dumper dumperFile] [-decision] [-angleList angleListMid] [-ec pathToFlat pathToPDB ecFile] [-ligand] [-angleBD angleBreakDownFile]
 
  Parameters:

    required:
  	PDBpathsFile 
	metal
	[shellCutoff shellElmt], elements are separated by capital letters. If not provided, cutoff is set based on atomic radius, and element includes everything but 'H'
	-i/d/dd/bs [(criteria threshold) statisticsFilePrefix], 
		-i: iterative process; -d: get bindinding ligand via chi-squared test; -dd: get binding ligand via single ligand test; -bs: just bootstrap
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
#         FILE:  "metalCoordinations.pl"
#
#        USAGE:  ./metalCoordinations.pl
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
use MPCGanalysis;

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

#use JSON;
#use Data::Dumper::Concise;

########## read in required arguments ##########
my $pathsFile = shift @ARGV;
my $metal = shift @ARGV;

my ($iaControl, $threshold, $shellCutoff, $shellElement, $shellsJsonFile); 
## bond length cutoff for the first shell in bootstrap
if ($ARGV[0] =~ /^\d+(\.\d+)?$/)
  {
  $shellCutoff = shift @ARGV;
  $shellElement = shift @ARGV;
  $shellsJsonFile = shift @ARGV if ($ARGV[0] !~ /\-/);
  }

unless ($ARGV[0] eq "-i" || $ARGV[0] eq "-d"|| $ARGV[0] eq "-bs" || $ARGV[0] eq "-dd")
  { print STDERR $help; exit;}
my $flow = shift @ARGV;

## This is not really -i/-d specific
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
elsif ($ARGV[0] =~ /^(n|non|nonModel)$/) ## nonModel is to remove the zinc shell for all later iterations if if didn't pass the probability threshold.
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

######################################## 
############# main process #############
########################################

my $analyzer = MPCGanalysis->new("pathsFile" => $pathsFile, 
				 "element" => uc($metal), 
				 "majorCGs" => ["Tetrahedral", "TrigonalBipyramidal", "Octahedral", "PentagonalBipyramidal"], 
				 "minLigNum" => 4,
				 "shellCutoff" => $shellCutoff,
				 "shellElement" => $shellElement,
				 "jsonFile" => $shellsJsonFile);
#print "$metal: ", $analyzer->{numCenter}, "\n";
#print "Cluster: ", $analyzer->{numCluster}, "\n";
#print "Usable: ", $analyzer->{usable}, "\n";
#print "Unusable: ", $analyzer->{unusable}, "\n\n";
#exit;

$analyzer->bootstrapCoordination($statOutFileName, $shellsJsonFile);

## required argument defines the workflow
if ($flow eq "-i") 
  { $analyzer->IAcoordination($statOutFileName, $iaControl, $threshold); }
elsif ($flow eq "-d") 
  { $analyzer->bindShellViaDist($statOutFileName); }
elsif ($flow eq "-dd")
  { $analyzer->shellViaAdjustDistStd($statOutFileName); }

########## Optional arguments ##########
my ($seqOpt, $seq, $header, $seqFile, $rInputOpt, $rInputFile, $rfOpt, $rfFile, $statsFile, $leaveOut, $jsonOpt, $jsonFile, $dumperOpt, $dumperFile, $decisionOpt, $angleListOpt, $angleListMid, $probInputFile);
my ($ECopt, $pathToFlat, $pathToPDB, $ecFile, $ligandOpt, $angleBreakDownOpt, $angleBreakDownFile, $bondLengthOpt, $bondLengthFile);
while (@ARGV) 
  {
  if ($ARGV[0] !~ /\-/)
    { print STDERR $help; exit;}
  my $option = shift @ARGV;

  if ($option eq "-s")
    {
    $seq = shift @ARGV;
    $header = shift @ARGV;
    $seqFile = shift @ARGV;
    
    unless ($seq eq "s" || $seq eq "a" || $seq eq "n")
      { print STDERR $help; exit;}
    unless ($header eq "b" || $header eq "s" || $header eq "ss" || $header eq "c")
      { print STDERR $help; exit;} 
    $analyzer->printSequences($seqFile, $seq, $header);
    }
  elsif ($option eq "-r")
    {
    $rInputFile = shift @ARGV;
    $leaveOut = shift @ARGV if ($ARGV[0] =~ /^(l|leave|leaveOut)$/);
    $statsFile = shift @ARGV if ($ARGV[0] !~ /\-/);
    &rPrint($analyzer, $rInputFile, $leaveOut, $statsFile) ;
    }
  elsif ($option eq "-rf")
    {
    $rfFile = shift @ARGV;
    $statsFile = shift @ARGV if ($ARGV[0] !~ /\-/);
    &rfPrint($analyzer, $rfFile, $statsFile) ;
    }
  elsif ($option eq "-probs")
    {
    $probInputFile = shift @ARGV; 
    $leaveOut = shift @ARGV;
    $statsFile = shift @ARGV if ($ARGV[0] !~ /\-/);
    &probPrint($analyzer, $probInputFile, $leaveOut, $statsFile) ;
    }
  elsif ($option eq "-bondLength")
    {
    $bondLengthFile = shift @ARGV;
    &blPrint($analyzer, $bondLengthFile);
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

###########################
## optional args' functions
###########################

sub rfPrint
  {
  my $analyzer = shift @_;
  my $rfFile = shift @_;
  my $statsFile = shift @_;

  open (RFID, ">", $rfFile) or die $!;
  my $stats = &readTableFile($statsFile) if ($statsFile);

  my %ids;
  map { my $ligNum = $_; map { $ids{$_->{shellObj}->metalID()} += 1; } (@{$analyzer->{coordinations}{$ligNum}}); } (keys %{$analyzer->{coordinations}});
  foreach my $ligNum (keys %{$analyzer->{coordinations}}) #("ten", "nine", "eight", "seven", "six", "five", "four")
    {
    foreach my $metalObj (@{$analyzer->{coordinations}{$ligNum}})
      {
      next if ($ids{$metalObj->{shellObj}->metalID()}) > 1;

      print RFID $metalObj->{shellObj}->metalID(), "\t";
      map { print RFID "$_\t";} ($metalObj->smallestAngle()); 
      print RFID $metalObj->{shellObj}->{center}->{method}, "\t";
      print RFID $metalObj->{shellObj}->{center}->{date}, "\t";
      print RFID $metalObj->{shellObj}->{center}->{resolution}, "\t";
      print RFID "\n";
      }
    }
  close RFID;
  }

sub rPrint
  {
  my $analyzer = shift @_;
  my $rInputFile = shift @_;
  my $leaveOut = shift @_;
  my $statsFile = shift @_;

  open (FID, ">", $rInputFile) or die $!;
  my $stats = &readTableFile($statsFile) if ($statsFile);
  my $stats = $analyzer->{stats} if ($flow eq "-i");

  #print FID "Metal_ID\tAngle_1\tAngle_2\tAngle_3\tAngle_4\tAngle_5\tAngle_6\tLigand_1\tLigand_2\tLigand_3\tLigand_4\t";
  #print FID "Bond_Length_1\tBond_length_2\tBond_length_3\tBond_length_4\tBi_status_1\tBi_status_2\tBi_status_3\tBi_status_4\tBi_status_5\tBi_status_6\t";
  #print FID "Method\tDate\tResolution\tB_factor_1\tB_factor_2\tB_factor_3\tB_factor_4\tOverall_bi_status\tLig_chain_ID_combo\tLig_residue_combo\tLig_atom_combo\tAmineN\t";
  #print FID "Chi_prob_Tet\tChi_prob_Bva\tChi_prob_Bvp\tChi_prob_Spv\tChi_prob_Spl\t" if ($statsFile);
  #print FID "\n";

  if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
    { $analyzer->calcChiCoordination($iaControl, $threshold, $leaveOut); }

  my %ids;
  map { my $ligNum = $_; map { $ids{$_->{shellObj}->metalID()} += 1; } (@{$analyzer->{coordinations}{$ligNum}}); } (keys %{$analyzer->{coordinations}});
  foreach my $ligNum (keys %{$analyzer->{coordinations}}) #("ten", "nine", "eight", "seven", "six", "five", "four")
    {
    foreach my $metalObj (@{$analyzer->{coordinations}{$ligNum}})
      {
      next if ($ids{$metalObj->{shellObj}->metalID()}) > 1;
      print FID $metalObj->{shellObj}->metalID(), "\t";

      print FID $metalObj->{shellObj}->{center}->{method}, "\t";
      print FID $metalObj->{shellObj}->{center}->{date}, "\t";
      print FID $metalObj->{shellObj}->{center}->{resolution}, "\t";

      print FID join (',', ($metalObj->orderedAngles())), "\t";
      map { print FID "$_\t";} ($metalObj->ligandAtomElement()); ## also have bond lengths for each ligand and bidentations for each angle
    
      print FID join (',', map {$_->{bFactor}} (@{$metalObj->{bestCombo}->{ligands}})), "\t";
      print FID $metalObj->bidentate(), "\t";
      map { print FID "$_\t";} ($metalObj->ligandCombos());

      #map { print FID $_, "\t";} (&coordProbs($metalObj->{shellObj}, $stats, $leaveOut)) if ($statsFile || $flow eq "-i");
      #map { print FID $_->{chiAngle}, "\t";} (@{$metalObj->{bestCombo}->{ligands}});
   
      print FID $metalObj->{shellObj}->{center}->{occupancy}, "\t"; 
      print FID $metalObj->{shellObj}->{center}->{solvent}, "\t";
      print FID "$ligNum\t", $metalObj->{bestCombo}->{probability}, "\t";
      
      print FID "\n";
      }
    }
  close FID;
  }

sub probPrint
  {   
  my $analyzer = shift @_;
  my $probInputFile = shift @_;
  my $leaveOut = shift @_;
  my $statsFile = shift @_;
  
  my $stats = ($statsFile)? &readTableFile($statsFile) : $analyzer->{stats};

  $analyzer->calcChiCoordination(0, 0, $leaveOut, $stats, $probInputFile);
  }



## Some other supplemental options
sub blPrint
  {
  my $analyzer = shift @_;
  my $bondLengthFile = shift @_;

  open (BLF, ">", $bondLengthFile) or die $!;
  print BLF join("\t", "metalID", "residueID", "element", "bondLength", "resolution", "rValue", "rFree"), "\n" ; 

  foreach my $cg (keys %{$analyzer->{coordinations}})
    {
    foreach my $metalObj (@{$analyzer->{coordinations}{$cg}})
      {
      my $comboLigands = $metalObj->{bestCombo}->{ligands};
      my $center = $metalObj->{shellObj}->{center};
      print BLF map {join("\t", $metalObj->{shellObj}->metalID(), $_->resID, $_->{element}, $center->distance($_), $_->{resolution}, $_->{rValue}, $_->{rFree}), "\n" ;} (@$comboLigands);
      }
    }
    close BLF;
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



if ($jsonOpt)
  {
  open (JOUT, '>', $jsonFile);
  my $jsonObj = JSON->new->allow_blessed->convert_blessed->encode( $analyzer->{nonModels} );
  print JOUT $jsonObj;
  }

if ($dumperOpt)
  {
  open (DOUT, '>', $dumperFile );
  print DOUT Dumper($analyzer->{nonModels});
  }

if ($decisionOpt)
  {
  foreach my $key ( sort keys %{$analyzer->{decisions}} )
    { print "$key, ", $analyzer->{decisions}{$key}, "\n"; }
  }

if ($angleListOpt)
  {
  foreach my $analyzerination (keys %{$analyzer->{coordinations}} )
    {
    next if (! $analyzer->{coordinations}{$analyzerination});

    my $fileName = "$analyzerination.$angleListMid.txt";
    open (FID, ">", $fileName) or die $!;
    foreach my $metalObj (@{$analyzer->{coordinations}{$analyzerination}})
      {
      map { print FID  "$_  \t";} ($metalObj->angleList());
      my $prob = $metalObj->{bestCombo}->{probability};
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
  foreach my $analyzerination (keys %{$analyzer->{coordinations}})
    {
    foreach my $metalObj (@{$analyzer->{coordinations}{$analyzerination}})
      {
      my $model = ref $metalObj;
      my $metalAtom = $metalObj->{shellObj}->{center};
      my $pdbid = $metalAtom->{PDBid};
      my $chainid = $metalAtom->{chainID};
      my $serial = $metalAtom->{residueNumber};

      my $id = "$pdbid.$chainid.$serial";
      my $inputId = lc($pdbid).uc($chainid);
  
      $$ecUnp{$id}{"uniprotID"} = [&getUnpidFromPdbid($inputId)];
      $$ecUnp{$id}{"coordination"} = $model;
      }  
    }

  while(my $input = <UNPFLAT>) #foreach uniprot id, check if it's metal-associated
    {
    chomp $input;
    next if ($input =~ /^\s*$/);
    next if ($input !~ /,/); #if the uniprot doesn't have ec # at all

    my %unpid = map {($_, 1);} (split (/; /, $input));
    foreach my $id (keys %$ecUnp)
      {
      if ( grep {$unpid{$_};} (@{$$ecUnp{$id}{"uniprotID"}}) ) # accession number matches => it's metal-associated
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
  foreach my $cg (values %{$analyzer->{coordinations}})
    {
    foreach my $metalObj (@$cg)
      {
      my @residues = map {$_->{residueName};} (@{$metalObj->{bestCombo}->{ligands}});
      $ligandCombo{join ".", sort {lc($a) cmp lc($b);} (@residues)}++;
      }
    }
  map {print "$_= ", $ligandCombo{$_}, "\n";} (sort {$ligandCombo{$b} <=> $ligandCombo{$a}} (keys %ligandCombo));
  }

if ($angleBreakDownOpt)
  {
  foreach my $cg (keys %{$analyzer->{rawAngles}})
    {
    foreach my $angel (keys %{$analyzer->{rawAngles}{$cg}} )
      {
      my $histogram=[];

      push @$histogram, scalar (grep {$_ > 0  && $_ < 15;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 15 && $_ < 25;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 25 && $_ < 35;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 35 && $_ < 45;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 45 && $_ < 55;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 55 && $_ < 65;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 65 && $_ < 75;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 75 && $_ < 85;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 85 && $_ < 95;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 95 && $_ < 105;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 105 && $_ < 115;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 115 && $_ < 125;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 125 && $_ < 135;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 135 && $_ < 145;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 145 && $_ < 155;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 155 && $_ < 165;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 165 && $_ < 175;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 175 && $_ < 180;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
      push @$histogram, scalar (grep {$_ > 180 ;} (@{$analyzer->{rawAngles}{$cg}{$angel}}) ) ;
 
      $analyzer->{rawAngles}{$cg}{$angel} = $histogram;
      }
    }

  &writeTableFile($angleBreakDownFile, $analyzer->{rawAngles});
  }


