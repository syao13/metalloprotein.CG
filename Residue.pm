## Residue.pm
#

package Residue;
use strict;
use Atom;

our @defaultDataMembers = (     "N" => 0,
				"CA" => 0,
                                "CB" => 0,
                                "CG" => 0
                                # "atom" => []
				);

sub new
  {
  my $class = shift @_;
  my $self = { @defaultDataMembers, "atoms" => [] };
  bless $self, ref $class || $class;

  $self->add(@_); # add an arbitrary number of atoms.
  return $self;
  }

sub add
  {
  my $self = shift @_;
  
  foreach my $atom (@_)
    {
    next if (@{$self->{atoms}} &&  $atom->resID() ne $self->residueID());

    push @{$self->{atoms}}, $atom; 
    $self->{N} = $atom if ($atom->{atomName} eq "N");
    $self->{CA} = $atom if ($atom->{atomName} eq "CA");
    $self->{CB} = $atom if ($atom->{atomName} eq "CB");
    $self->{CG} = $atom if ($atom->{atomName} eq "CG" || $atom->{atomName} eq "SG" || $atom->{atomName} eq "OG" || $atom->{atomName} eq "OG1" || $atom->{atomName} eq "CG1");
    }
  }
    

sub residueID
  {
  my $self = shift @_;

  return $self->{atoms}[0]->resID();
  }

sub cmp
  {
  my $self = shift @_;
  my $other = shift @_;

  my $selfAtom = $self->{atoms}[0];
  my $otherAtom = $other->{atoms}[0];
  
  return ($selfAtom->{chainID} eq $otherAtom->{chainID}) ? ($selfAtom->{residueNumber} <=> $otherAtom->{residueNumber}) : ($selfAtom->{chainID} cmp $otherAtom->{chainID});
  }

sub chiOneAngle
  {
  my $self = shift @_;

  return 0 if (! ($self->{N} && $self->{CA} && $self->{CB} && $self->{CG}));

  return $self->{N}->torsionAngle($self->{CA}, $self->{CB}, $self->{CG});
  }

