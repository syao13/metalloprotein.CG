## Sequence.pm
###
package Sequence;
use strict;

## Constructor: new
### Methods: 

our @defaultSequence = (
                "header" => "",
		"sequence" => ""
                 );


## Convert arrays of 3-letter aa sequences from SEQRES section into one 1-letter sequence
sub createSEQRESseq
  {
  my $class = shift @_;  
  $class = ref $class || $class;

  my $chain = shift @_;
  my $seqs = shift @_;
  my $atoms = shift @_;

  foreach my $atom (@$atoms)
    {
    next if ($atom->{recordType} ne "ATOM");

    if ( $chain eq $atom->{chainID} )
      { 
      $chain = $chain."|".$chain.$atom->{residueNumber};
      last;
      }
    }

  my $seq1letter;
  foreach my $seq3letter (@$seqs)
    {
    my @aa3letter = split(/ /, $seq3letter);

    $seq1letter = $seq1letter.join('', (map {(&_aaCode($_))? &_aaCode($_): "X" ;} (@aa3letter))); 
    }

  return $class->new("header" => $chain, "sequence" => $seq1letter);
  }


## Generate sequence from ATOM object array
sub createATOMseqs
  {
  my $class = shift @_;
  $class = ref $class || $class;

  my $atoms = shift @_;
  my %chainSeqs;
  my %idRep; ## use each chain_ID.residue_number only once
  my %chainFirstNum; ## get the first residue number for each chain
  my %chainNum; ## count what residue number each chain get to, the missing ones will be denote as "U";
  foreach my $atom (@$atoms)
    {
    next if ($atom->{recordType} ne "ATOM");

    # acquire the first number for each chain
    my $chain = $atom->{chainID};
    $chainFirstNum{$chain} = $chain.$atom->{residueNumber} if (! $chainFirstNum{$chain});
    $chain = $chain."|".$chainFirstNum{$chain};

    # use each residue only once
    my $id = $atom->resID();
    next if $idRep{$id} == 3;

    # insert "U" to represent missing residue
    while ($chainNum{$chain} = ($chainNum{$chain})? ($chainNum{$chain}+1) : $atom->{residueNumber})
      {
      if ($atom->{residueNumber} == $chainNum{$chain})
	{last;}
      elsif ($atom->{residueNumber} < $chainNum{$chain})
	{
	$chainNum{$chain} = $atom->{residueNumber}; 
	$chainSeqs{$chain} = substr($chainSeqs{$chain}, 0, ($atom->{residueNumber} - substr($chain, 3)));
	
	last;
	}
      else
        { 
	$chainSeqs{$chain} = $chainSeqs{$chain}."-"; 
	if ($chainNum{$chain} > 9999)
	  {
	  print $atom->{PDBid}, ": flag, ATOM sequence error!\n"; 
	  last;
	  }
	}
      }

    my $oneLcode = (&_aaCode($atom->{residueName}))? &_aaCode($atom->{residueName}) : "X";
    $chainSeqs{$chain} = $chainSeqs{$chain}.$oneLcode;
    $idRep{$id} = 3;
    }

  map {$_ =~ s/X+$//; } (values %chainSeqs); # remove Xs at the end from HETATMs
  return [ map { $class->new("header" => $_, "sequence" => $chainSeqs{$_}); } (keys %chainSeqs) ];
  }


## Generate sequence from ATOM object array
sub createATOMnum
  {
  my $class = shift @_;
  $class = ref $class || $class;

  my $atoms = shift @_;
  my %chainSeqs;
  my %idRep; ## use each chain_ID.residue_number only once
  my %chainFirstNum; ## get the first residue number for each chain
  my %chainNum; ## count what residue number each chain get to, the missing ones will be denote as "U";

  foreach my $atom (@$atoms)
    {
    next if ($atom->{recordType} ne "ATOM");

    # acquire the first number for each chain
    my $chain = $atom->{chainID};
    $chainFirstNum{$chain} = $chain.$atom->{residueNumber} if (! $chainFirstNum{$chain});
    $chain = $chain."|".$chainFirstNum{$chain};

    # use each residue only once
    my $id = $atom->resID();
    next if $idRep{$id} == 3;

    # insert "U" to represent missing residue    
    while ($chainNum{$chain} = ($chainNum{$chain})? ($chainNum{$chain}+1) : $atom->{residueNumber})
      {
      if ($atom->{residueNumber} == $chainNum{$chain})
        {last;}
      elsif ($atom->{residueNumber} < $chainNum{$chain})
        {
        $chainNum{$chain} = $atom->{residueNumber};
	my @numbers = split (/\./, $chainSeqs{$chain});
	my @toCurrNum = splice (@numbers, 0, ($atom->{residueNumber} - substr($chain, 3)));
	$chainSeqs{$chain} = join (".", @toCurrNum);

        last;
        }
      else
        {
        $chainSeqs{$chain} = $chainSeqs{$chain}.".".$chainNum{$chain};
        if ($chainNum{$chain} > 9999)
          {
          print $atom->{PDBid}, "flag\n";
          last;
          }
        }
      }

    $chainSeqs{$chain} = ($chainSeqs{$chain})? ($chainSeqs{$chain}.".".$chainNum{$chain}) : ($chainNum{$chain});
    $idRep{$id} = 3;
    }

  map {$_ =~ s/X+$//; } (values %chainSeqs); # remove Xs at the end from HETATMs
  return [ map { $class->new("header" => $_, "sequence" => $chainSeqs{$_}); } (keys %chainSeqs) ];
  }
  
    
 
sub new
  {
  my $class = shift @_;
  my $self = { @defaultSequence, @_};

  return bless $self, ref $class || $class;
  }


sub updateHeader
  {
  my $self = shift @_;

  my $firstAtomNum = (split ('\|', $self->{header}))[-1];

  $self->{header} = shift @_;  
  while (my $item = shift @_)
    {$self->{header} = $self->{header}."|".$item; }
  $self->{header} = $self->{header}."|".$firstAtomNum;

  }


sub printFasta
  {
  my $self = shift @_;
  my $fh = shift @_;  ## file handle
  my $len = (@_)? shift @_: 80;

  my $seq = $self->{sequence};
  print $fh ">".$self->{header}, "\n" ;
  while (my $chunk = substr($seq, 0, $len, "")) 
    { print $fh "$chunk\n"; }
  }


sub printNum
  {
  my $self = shift @_;
  my $fh = shift @_;  ## file handle

  print $fh ">".$self->{header}, "\n" ;
  print $fh $self->{sequence}, "\n"; 
  }



sub _aaCode
  {
  my $aa3 = shift @_;
  
  my %code = (
	"ALA"     =>   "A",
        "ARG"	  =>   "R",
	"ASN"	  =>   "N",
	"ASP"	  =>   "D",
	"CYS"     =>   "C", 
        "GLU"     =>   "E",
        "GLN"     =>   "Q",
        "GLY"     =>   "G",
        "HIS"     =>   "H",
        "ILE"     =>   "I",
        "LEU"     =>   "L",
        "LYS"     =>   "K",
        "MET"     =>   "M",
        "PHE"     =>   "F",
        "PRO"     =>   "P",
        "SER"     =>   "S",
        "THR"     =>   "T",
        "TRP"     =>   "W",
        "TYR"     =>   "Y",
        "VAL"     =>   "V",
        "ASX"     =>   "B",
        "GLX"     =>   "Z",  
        "..."     =>   "X");

  return $code{$aa3};
  }     



