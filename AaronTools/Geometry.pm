#Contributors Yanfei Guan and Steven E. Wheeler

use lib $ENV{'PERL_LIB'};
use lib $ENV{'AARON'};

use Constants qw(:PHYSICAL :COMPARE);
use AaronTools::Atoms qw(:BASIC :LJ);
use AaronTools::Molecules;
use AaronTools::FileReader;

our $boltzman = BOLTZMANN;
our $CONNECTIVITY = CONNECTIVITY;
our $UNROTATABLE_BOND = UNROTATABLE_BOND;
our $CUTOFF = CUTOFF;
our $mass = MASS;
our $radii = RADII;
our $AARON = $ENV{'AARON'};
our $rij = RIJ;
our $eij = EIJ;

package AaronTools::Geometry;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
use Math::MatrixReal;
use Data::Dumper;

sub new {
    my $class = shift;
    my %param = @_;
    my $self = {
        name => $param{name},
        elements => $param{elements},
        flags => $param{flags},
        coords => $param{coords},
        connection => $param{connection},
        constraints => $param{constraints},
    };

    bless $self, $class;

    $self->{elements} //= [];
    $self->{flags} //= [];
    $self->{coords} //= [];
    $self->{connection} //= [];
    $self->{constraints} //= [];

    unless (@{$self->{coords}}) {
        if ($self->{name}) {
            if (-f "$self->{name}.xyz") {
                $self->read_geometry("$self->{name}.xyz");
            }elsif ($self->{name} && AaronTools::Molecules::built_in($self->{name})) {
                $self->{coords} = AaronTools::Molecules::coords($self->{name});
                $self->{elements} = AaronTools::Molecules::elements($self->{name});
                $self->{flags} = AaronTools::Molecules::flags($self->{name});
            }
        }
    }


    $self->refresh_connected() unless @{$self->{connection}};

    return $self;
}


sub set_name {
    my ($self, $name) = @_;
    $self->{name} = $name;
}


sub name {
    my $self = shift;
    return $self->{name};
}


sub elements {
    my $self = shift;
    return $self->{elements};
}


sub flags {
    my $self = shift;
    return $self->{flags};
}


sub coords {
    my $self = shift;
    return Math::MatrixReal->new_from_rows($self->{coords});
}


sub connection {
    my $self = shift;
    return $self->{connection};
}


sub constraints {
    my $self = shift;
    return $self->{constraints};
}


sub set_constraints {
    my ($self, $constraints) = @_;
    $self->{constraints} = $constraints;
}


sub freeze_atoms {
    my ($self, $atoms) = @_;
    for my $atom (@$atoms) {$self->{flags}->[$atom] = -1;};
}


sub copy {
    my $self = shift;
    my $new =  new AaronTools::Geometry( name => $self->{name},
                                         elements => [ @{ $self->{elements} } ],
                                         flags => [ @{ $self->{flags} }],
                                         coords => [ map { [ @$_ ] } @{ $self->{coords} } ],
                                         connection => [ map { [ @$_ ] } @{ $self->{connection} } ],
                                         constraints => [ map { [ @$_ ] } @{ $self->{constraints} } ] );
};


sub update_coords {
    my $self = shift;
    my %param = @_;
    
    for my $i (0..$#{ $param{targets} }) {
        $self->{coords}->[$param{targets}->[$i]] = $param{coords}->[$i];
    }
}


sub append {
    my ($self, $geo_2) = @_;
    my $num_self = $#{ $self->{elements} } + 1;

    $self->{elements} = [@{ $self->{elements} }, @{ $geo_2->{elements} }];
    $self->{coords} = [@{ $self->{coords} }, @{ $geo_2->{coords} }];
    $self->{flags} = [@{ $self->{flags} }, @{ $geo_2->{flags} }];

    my @connection_2_new = map{ [map {$_ + $num_self} @$_] }  @{ $geo_2->{connection} };
    $self->{connection} = [@{ $self->{connection} }, @connection_2_new];
}


sub subgeo {
    my ($self, $groups) = @_;
    $groups //= [0..$#{ $self->{elements} }];

    my @elements;
    my @flags;
    my @coords;
    my @connection;

    for my $atom (@$groups) {
        push (@elements, $self->{elements}->[$atom]);
        push (@flags, $self->{flags}->[$atom]);
        push (@coords, $self->{coords}->[$atom]);
    }
    my $subgeo = new AaronTools::Geometry( name => $self->{name},
                                           elements => [@elements],
                                           flags => [@flags],
                                           coords => [@coords]);

    return $subgeo;
}


sub splice_atom {
    my ($self, $target, $number_splice, $geo_2) = @_;
    
    for my $entry ('elements', 'flags', 'coords', 'connection') {
        if ($geo_2) {
            splice(@{ $self->{$entry} }, $target, $number_splice, @{ $geo_2->{$entry} });
        }else {
            splice(@{ $self->{$entry} }, $target, $number_splice);
        }
    };
}


sub delete_atom {
    my ($self, $groups) = @_;
    $groups //= [];

    for my $atom (sort { $b <=> $a } @$groups) {
        $self->splice_atom($atom, 1);
    }
}


#Removes fragment by cutting $atom1-$atom2 bond and removing $atom2 and all connected atoms, replaces with hydrogen along original $atom1-$atom2 bond.
#remove_fragment(\@coords, $atom1, $atom2);
sub remove_fragment {
    my ($self, $atom1, $atom2) = @_;

    my $fragment = $self->get_all_connected($atom2, $atom1);
    my @fragment_no_atom2 = grep { $_ != $atom2 } @$fragment;
    $self->{elements}->[$atom2] = 'H';
    $self->correct_bond_length( atom1 => $atom1, atom2 => $atom2 );
    $self->delete_atom([@fragment_no_atom2]);
    $self->get_connected();
}


sub clean_structure {
    my ($self) = @_;
    my @phantom;
    for my $i (0..$#{ $self->{elements} }) {
        if ($self->{elements}->[$i] eq 'X') {
            push (@phantom, $i);
        }
    }
    $self->delete_atom([@phantom]);
    $self->refresh_connected();
}


#makes a copy of a @coords array, mirroring x-values.  Returns new @coords array
#mirror_coords($ref_to_coords);
sub mirror_coords {
    my ($self, $axis) = @_;

    $axis //= '';

    SWITCH: {
        if ($axis =~ /[yY]/) {
            foreach my $atom (0..$#{ $self->{elements} }) {
                $self->{coords}->[$atom]->[1] *= -1;
            }
            last SWITCH;
        }
        if ($axis =~ /[zZ]/) {
            foreach my $atom (0..$#{ $self->{elements} }) {
                $self->{coords}->[$atom]->[2] *= -1;
            }
            last SWITCH;
        }
        foreach my $atom (0..$#{ $self->{elements} }) {
            $self->{coords}->[$atom]->[0] *= -1;
        }
    }
}


sub read_geometry {
    my $self = shift;
    my ($file) = @_;
    my ($elements, $flags, $coords, $constraints) = AaronTools::FileReader::grab_coords($file);

    $self->{elements} = $elements;
    $self->{flags} = $flags;
    $self->{coords} = $coords;

    if ($constraints) {
        my @constraints;
        for my $constraint (@$constraints) {
            my $d = $self->distance( atom1 => $constraint->[0],
                                     atom2 => $constraint->[1] );
            push(@constraints, [$constraint, $d]);
            $self->{constraints} = [@constraints];
        }
    }

    $self->refresh_connected();
    return 1;
}


#get connected for subgroups atoms. If subgroups are not provided,
#get connected for all atoms,
#The subgroups are used when make substituents.
sub get_connected {
  my ($self, $subgroups) = @_;
  $subgroups //= [0..$#{ $self->{elements} }];

  my $tolerance = 0.2;

  my $connected = [];
  foreach my $atom1 (@$subgroups) {
    my @row;
    foreach my $atom2 (@$subgroups) {
      my $distance = $self->distance(atom1 => $atom1, atom2 => $atom2);

      my $cutoff = $radii->{$self->{elements}->[$atom1]} + 
                   $radii->{$self->{elements}->[$atom2]} + 
                   $tolerance;
      if($distance < $cutoff && $atom1 != $atom2) {
        push(@row, $atom2);
      }
    }
    $self->{connection}->[$atom1] = union($self->{connection}->[$atom1], [@row]);
  }
} #end of get_connected


sub refresh_connected {
  my $self = shift;
  my $subgroups = [0..$#{ $self->{elements} }];

  my $tolerance = 0.2;

  my $connected = [];
  foreach my $atom1 (@$subgroups) {
    my @row;
    foreach my $atom2 (@$subgroups) {
      my $distance = $self->distance(atom1 => $atom1, atom2 => $atom2);

      my $cutoff = $radii->{$self->{elements}->[$atom1]} + 
                   $radii->{$self->{elements}->[$atom2]} + 
                   $tolerance;
      if($distance < $cutoff && $atom1 != $atom2) {
        push(@row, $atom2);
      }
    }
    $self->{connection}->[$atom1] = [@row];
  }
} #end of get_connected



sub get_all_connected {
  my ($self, $start_atom, $avoid_atom) = @_;
  my @connected = @{$self->{connection}};
  my @connected_temp = map { [@$_] } @connected;
  #remove avoid_atom from list of atoms connected to $start_atom
  my $avoid_atom_found = 0;
  foreach my $atom (0..$#{$connected_temp[$start_atom]}) {
    if($connected_temp[$start_atom][$atom] == $avoid_atom) {
      $connected_temp[$start_atom][$atom] = -1;
      $avoid_atom_found = 1;
    }
  }

  my @positions = ($start_atom);
  #I can probably use something simpler than a hash here (since I'm not really using the values)
  my %visited = ( $start_atom => '0') ;	#keys are numbers of the visited atoms, values are not used

  #loop until @positions is empty
  while(@positions) {
    my $position = shift(@positions);
    foreach my $atom (@{$connected_temp[$position]}) { 	#grab all atoms connected to current atom and add to queue (unless already visited)
      if($atom >= 0 && !exists $visited{$atom}) {
        push(@positions, $atom);
        $visited{$atom} = 0;
      }
    }
  }

  my @all_connected_atoms = keys %visited;
  #Check all_connected_atoms for avoid_atom
  foreach (@all_connected_atoms) {
    if($_ == $avoid_atom) {
      return ();
    }
  }
  #change the start_atom in the @all_connected_atoms to the first element
  if ($all_connected_atoms[1]) {
    @all_connected_atoms = grep {$_ != $start_atom} @all_connected_atoms;
    @all_connected_atoms = sort {$a <=> $b} @all_connected_atoms;
    unshift @all_connected_atoms, $start_atom;
  }
  return [@all_connected_atoms];
} #End sub get_all_connected


sub check_connectivity {
    my ($self) = @_;

    for my $atom (0..$#{ $self->{elements} }) {
        if ($#{ $self->{connection}->[$atom] } + 1 > 
            $CONNECTIVITY->{$self->{elements}->[$atom]}) {return $atom;}
    }
    return -1;
}


#put in the start and end atom, return what a substituent object.
#if the substituent is built-in, information about conformation will be
#included. Otherwise, only geometries.
sub detect_substituent {
    my $self = shift;
    my %params = @_;

    my ($target, $end) = ($params{target}, $params{end});

    my $sub_atoms = $self->get_all_connected($target, $end);

    my $substituent = $self->subgeo($sub_atoms);
    $substituent->refresh_connected();

    bless $substituent, "AaronTools::Substituent";

    $substituent->compare_lib();
    $substituent->{end} = $end;
    return $substituent;
}


sub examine_constraints {
    my $self = shift;

    my $return = 0;
    my $con_return;
    for my $constraint (@{ $self->{constraints} }) {
        my $bond = $constraint->[0];
        my $d_con = $constraint->[1];

        my $d = $self->distance( atom1 => $bond->[0],
                                 atom2 => $bond->[1] );
        if ($d - $d_con > $CUTOFF->{D_CUTOFF}) {
            $return = 1;
            $con_return = $bond;
            last;
        }elsif ($d_con - $d > $CUTOFF->{D_CUTOFF}) {
            $return = -1;
            $con_return = $bond;
            last;
        }
    }

    return ($return, $con_return);
}






sub distance {
    my $self = shift;
    my %param = @_;
    my ($atom1, $atom2, $geometry_2) = ( $param{atom1},
                                         $param{atom2},
                                         $param{geometry2} );

    my $bond = $self->get_bond($atom1, $atom2, $geometry_2);
    my $distance = abs($bond);
    return $distance;
}


sub change_distance {
    my $self = shift;
    my %params = @_;

    my ($atom1, $atom2, $distance, 
        $by_distance,
        $move_atom2, $move_frag) = ( $params{atom1},
                                     $params{atom2},
                                     $params{distance},
                                     $params{by_distance},
                                     $params{fix_atom1},
                                     $params{translate_group} );

    $move_atom2 //= 0;
    $move_frag //= 1;

    my ($all_connected_atoms1, $all_connected_atoms2);
    if ($move_frag) {
        $all_connected_atoms2 = $self->get_all_connected($atom2, $atom1);
        if ($move_atom2) {
            $all_connected_atoms1 = $self->get_all_connected($atom1, $atom2);
        }
    }else {
        $all_connected_atoms2 = [$atom2];
        if ($move_atom2) {
            $all_connected_atoms1 = [$atom1];
        }
    }

    $self->_change_distance( all_atoms1 => $all_connected_atoms1,
                             all_atoms2 => $all_connected_atoms2,
                             atom1 => $atom1,
                             atom2 => $atom2, 
                             distance => $distance,
                             by_distance => $by_distance );
}


sub _change_distance {
    my $self = shift;
    my %params = @_;

    my ($all_atoms1, $all_atoms2, 
        $atom1, $atom2, $distance,
        $by_distance) = ($params{all_atoms1}, $params{all_atoms2},
                         $params{atom1}, $params{atom2}, $params{distance},
                         $params{by_distance});

    my $current_distance = $self->distance(atom1 => $atom1, atom2 => $atom2);
    my $difference;
    if ($distance) {
        $difference = $distance - $current_distance;
    }elsif($by_distance) {
        $difference = $by_distance;
    }

    my $v12 = $self->get_bond($atom2, $atom1);

    my ($v1, $v2);
    if($all_atoms1) {
        $v2 = $v12 * $difference / $current_distance;
        $v1 = V(0, 0, 0);
    }else {
        $v2 = $v12 * $difference / (2*$current_distance);
        $v1 = -$v12 * $difference / (2*$current_distance);
    }

    $self->coord_shift($v1, $all_atoms1);
    $self->coord_shift($v2, $all_atoms2);
}


sub correct_bond_length {
    my $self = shift;
    my %params = @_;
    my ($atom1, $atom2) = ($params{atom1}, $params{atom2});
    my $new_distance = $radii->{$self->{elements}->[$atom1]} 
                       + $radii->{$self->{elements}->[$atom2]};

    $self->change_distance( atom1 => $atom1,
                            atom2 => $atom2,
                            distance => $new_distance,
                            fix_atom1 => 1,
                            translate_group => 1 );
}


sub coord_shift {
    my ($self, $v, $targets) = @_;
    $targets //= [0..$#{ $self->{coords} }];

    foreach my $atom (@$targets) {
        my $new = $self->{coords}->[$atom] + $v;
        #This is to maintain the array ref
        for my $i (0..$#{ $new }) {$self->{coords}->[$atom]->[$i] = $new->[$i]};
    }
}


sub quat_rot {
    my ($self, $a, $w, $targets) = @_;

    $targets //= [ 0..$#{ $self->{elements} }];

    for my $atom (@$targets) {
        my $vec = $self->get_point($atom);
        my $wx = $w x $vec;
        my $new_vec = $vec + 2*$a*$wx + 2*($w x $wx);
        #This is to maintain the array ref
        for my $i (0..$#{ $new_vec }) {$self->{coords}->[$atom]->[$i] = $new_vec->[$i]};
    }
}


sub genrotate {
    my ($self, $v, $angle, $targets) = @_;

    my $a = cos($angle/2);
    $v /= $v->norm();
    $v *= sin($angle/2);

    $self->quat_rot($a, $v, $targets);
}


sub rotate {
    my ($self, $axis, $angle) = @_;
    Switch:{
        if ($axis =~ /[Xx]/) {
            $self->genrotate(V(1,0,0), $angle);
            last Switch;
        }
        if ($axis =~ /[Yy]/) {
            $self->genrotate(V(0,1,0), $angle);
            last Switch;
        }
        if ($axis =~ /[Zz]/) {
            $self->genrotate(V(0,0,1), $angle);
            last Switch;
        } 
        die "Can only rotate around Cartesian axes!\n";
    }
}


sub center_genrotate {
    my ($self, $point, $v, $angle, $targets) = @_;
    
    my $shift_v = $point =~ /^\d+$/ ?
                  $self->get_point($point) : $point;

    $self->coord_shift($shift_v*-1, $targets);
    $self->genrotate($v, $angle, $targets);
    $self->coord_shift($shift_v, $targets);
}


sub angle {
    my ($self, $atom1, $atom2, $atom3) = @_;
    
    my $bond1 = $self->get_bond($atom1, $atom2);
    my $bond2 = $self->get_bond($atom3, $atom2);

    my $angle = atan2($bond1, $bond2);

    return $angle;
}


sub dihedral {
    my($self, $atom1, $atom2, $atom3, $atom4) = @_;

    my $bond12 = $self->get_bond($atom1, $atom2);
    my $bond23 = $self->get_bond($atom2, $atom3);
    my $bond34 = $self->get_bond($atom3, $atom4);

    my $dihedral = atan2((($bond12 x $bond23) x ($bond23 x $bond34)) * $bond23 /
                          abs($bond23), ($bond12 x $bond23) * ($bond23 x $bond34));

    return rad2deg($dihedral);
}


#changes dihedral about bond atom1-atom2 by an angle angle_change (in radians!)
#change_dihedral(atom1, atom2, angle_change, ref_to_coords
sub change_dihedral {
    my ($self, $atom1, $atom2, $angle) = @_;
    
    my $connected_atoms = $self->get_all_connected($atom1, $atom2);

    if ($#{ $connected_atoms } < 0) {
        print "Cannot change this dihedral angle...\n";
        return 0;
    }

    my $bond = $self->get_bond($atom2, $atom1);

    $self->center_genrotate($atom1, $bond, $angle, $connected_atoms);
}


sub set_dihedral {
    my ($self, $atom1, $atom2, $atom3, $atom4, $new_tau) = @_;

    my $tau = $self->dihedral($atom1, $atom2, $atom3, $atom4);
    $self->change_dihedral($atom2, $atom3, deg2rad($new_tau - $tau));
}


#centers molecule so that average of atom1, atom2, and atom3 is at origin, and orients molecule so that atom1 is on x-axis and atom2 and atom3 are as close as possible to XY-plane
#center_ring($atom1, $atom2, $atom3, $ref_to_coords);
sub center_ring {
    my ($self, $atom1, $atom2, $atom3) = @_;

    my $com = V(0, 0, 0);
    for my $atom ($atom1, $atom2, $atom3) {
        $com += $self->get_point($atom);
    }
    $com /= 3; 
    
    #shift geom to COM
    $self->coord_shift(-1*$com);
    
    #Put atom1 along x-axis
    my $v1 = $self->get_point($atom1);
    my $vx = V(1,0,0);
    my $cross1 = $v1 x $vx;
    my $angle1 = atan2($v1, $vx);
    $self->genrotate($cross1, -$angle1);
    
    #Now put atom2 and atom3 in XY-plane (or as close as possible)
    my ($v2, $v3) = ($self->get_point($atom2), $self->get_poing($atom3));
    $v2->[0] = 0;
    $v3->[0] = 0;
    my $vz = V(0,0,1);
    my $cross2 = $v2 x $vz;
    my $cross3 = $v3 x $vz;
    my $angle2;
    if($v2->norm() != 0) {
      $angle2 = asin($cross2->norm()/$v2->norm());
    } else {
      $angle2 = pi()/2;
    }
    my $angle3;
    if($v3->norm() != 0) {
      $angle3 = asin($cross3->norm()/$v3->norm());
    } else {
      $angle3 = pi()/2;
    }
    my $av_angle = pi()/2 - ($angle2*2)/2;

    $self->genrotate(V(1, 0, 0), -$av_angle);
}


#rotate substitute
#sub_rotate(target,angle,\@coords)
sub sub_rotate{
   my ($self, $target, $angle) = @_;

   my ($atom2, $targets) = $self->get_sub($target);
   my $axis = $self->get_bond($target, $atom2);
   $self->center_genrotate($target, $axis, $angle, $targets);
} #End sub sub_rotate


sub get_point {
    my ($self, $atom) = @_;
    my $vector = V(@{ $self->{coords}->[$atom] });
    return $vector;
}

#This is to find the substituent for a given atom
sub get_sub{
    my ($self, $target) = @_;

    my @nearst = @{$self->{connection}->[$target]};
    #now we need to define which nearst atom is the target.
    my $min = 999;
    my $atom2;
    my $targets;
    for (@nearst) {
        my $get_all = $self->get_all_connected($target, $_);
        if ($get_all && $#{ $get_all } < $min) {
          $min = $#{ $get_all };
          $atom2 = $_;
          $targets = $get_all;
        }
    }
    return($atom2, $targets);
}


#Replaces atom $target with substituent $sub (name of XYZ file)
#Returns coords in same orientation as given, with atom ordering preserved
#&substitute($target, $sub, \@coords, $no_min);
#TO DO: after adding substituent, rotation added atoms around new bond to maxmimize distance to all other atoms!
sub substitute {
    my $self = shift;
    my %param = @_;

    my ($target, $sub, $minimize_torsion) = ( $param{target},
                                              $param{sub},
                                              $param{minimize_torsion} );

    my ($end, $old_sub_atoms) = $self->get_sub($target);  

    my $sub_object = new AaronTools::Substituent( name => $sub, end => $end );

    $self->_substitute( old_sub_atoms => $old_sub_atoms,
                                  sub => $sub_object,
                                  end => $end,
                     minimize_torsion => $param{minimize_torsion} );
} #End sub substitute


#This is to substitute the atoms directly by providing a set of atoms number
#The constraint and to_sub will be modiyed.
#Don't call this unless you new what you are doing.
#For general users use the substitute instead by providing target and sub name.
sub _substitute {
    my $self = shift;
    my %params = @_;

    my ($old_sub_atoms, $sub, $end, $minimize_torsion) = ( $params{old_sub_atoms},
                                                           $params{sub},
                                                           $params{end},
                                                           $params{minimize_torsion} );
    my $target = $old_sub_atoms->[0];

    $sub->_align_on_geometry( geo => $self,
                           target => $target,
                              end => $end );

    
    #replace target with first atom of substituent
    $self->splice_atom($target, 1, $sub->subgeo([0])->copy());

    #if any element remaining in the @targets, remove them
    my $delete_atoms=[];
    for my $i(1..$#{ $old_sub_atoms }) {$delete_atoms->[$i-1] = $old_sub_atoms->[$i]};
    $target -= grep { $_ < $target } @$delete_atoms;
    $end -= grep { $_ < $end } @$delete_atoms;

    #modify the constraint, since the deleted atoms can change atom numbers
    $self->_rearrange_con_sub( delete_atoms => $delete_atoms);

    $self->delete_atom($delete_atoms);
    $self->refresh_connected();
    
    #build list of substituent coords
    my $old_num_atoms = $#{ $self->{elements} } + 1;
    $self->append($sub->subgeo([1..$#{ $sub->{elements} }])->copy());

    #modify the distance between end and target
    $self->get_connected([$target, ($old_num_atoms..$#{ $self->{elements} })]);
    $self->correct_bond_length( atom1 => $end, atom2 => $target );

    if ($minimize_torsion) {
        $self->minimize_torsion(start_atom => $target, 
                                  end_atom => $end);
    }
}


sub fused_ring {
  my ($self, $target1, $target2, $type) = @_;

  #get connected atoms
  my @connected = @{ $self->{connection} };
  if($#{$self->{connection}->[$target1]} > 0 
     || $#{$self->{connection}->[$target2]} > 0) {
    print "Trying to substitute non-monovalent atom!\n";
    return 0;
  } 

  #Figure out what needs to be added to match fused ring
  my $path_length = $self->shortest_path($target1, $target2);
  if($path_length < 3 || $path_length > 5) {
    print "Can't figure out how to build fused ring connecting those two atoms...\n";
    return 0;
  }
  #check to make sure known type
  if($type !~ /a_pinene/ && $type !~ /LD_chair/ && $type !~ /LU_chair/ 
     && $type !~ /D_boat/ && $type !~ /U_boat/ && $type !~ /Ar/) {
    print "Unknown ring type!\n";
    return 0;
  }

  my $ring;

  if ($type =~ /Ar/) {
    $ring = new AaronTools::Geometry( name => 'Ar_ring' );
    $ring->read_geometry("$AARON/Ring_fragments/six_$path_length.xyz");
  }
  else {	#Chairs
    if($path_length != 3) {
      print "Can't figure out how to build fused ring connecting those two atoms...\n";
      return 0;
    } else {
      $ring = new AaronTools::Geometry( name => 'Chairs' );
      $ring->read_geometry("$AARON/Ring_fragments/Chairs/$type.xyz");
    }
  }
  
  #get nearest neighbors from @connected
  my $nearest_neighbor1 = $self->{connection}->[$target1]->[0];
  my $nearest_neighbor2 = $self->{connection}->[$target2]->[0];

  #shift coords so that nearest_neighbor1 is at the origin
  my $origin = $self->get_point($nearest_neighbor1);
  $self->coord_shift(-1*$origin);

  #Orient so that nearest_neighbor2 is along x-axis
  my $v1 = $self->get_point($nearest_neighbor1);
  my $v2 = $self->get_point($nearest_neighbor2);

  my $vx = V(1,0,0);
  my $vz = V(0,0,1);
  my $cross1 = $v2 x $vx;
  my $angle1 = atan2($v2, $vx);
   #accounts for sign of dot product to get sign of angle right!
  $self->genrotate($cross1, $angle1) unless (abs($cross1) == 0);

  #final rotation around x-axis to bring $target1 and $target2 into XY-plane
  $v1 = $self->get_point($target1);
  $v2 = $self->get_point($target2);
  my $chi1 = atan2($v1->[2], $v1->[1]);
  my $chi2 = atan2($v2->[2], $v2->[1]);
  my $chi = ($chi1 + $chi2)/2;
  $self->rotate('x',-$chi);

  my @sub_atoms = ($target1, $target2);
  my $coord_num_old = $#{ $self->{elements} };
  
  #replace target1 with 1st ring atom
  #replace target2 with 2nd ring atom
  #Add remainder of ring coords
  my $ring1 = $ring->subgeo([0])->copy();
  my $ring2 = $ring->subgeo([1])->copy();
  my $ring_remain = $ring->subgeo([2..$#{ $ring->{elements} }])->copy();
  $self->splice_atom($target1, 1, $ring1);
  $self->splice_atom($target2, 1, $ring2);
  $self->append($ring_remain->copy());
  push(@sub_atoms, $coord_num_old+1..$#{ $self->{elements} });
  
  #Return geometry to original orientation/position
  $self->rotate('x',$chi);
  $self->genrotate($cross1, -$angle1) unless (abs($cross1) == 0);
  $self->coord_shift($origin);

  return [@sub_atoms];
} #End sub fused_ring

#calculates  LJ-6-12 potential energy based on autodock Rij and Eij parameters
#simply ignores any atom pair involving elements for which parameters are missing (which shouldn't be anything!).
sub LJ_energy {
  my ($self) = @_;

  my $energy = 0;

  foreach my $atom1 (0..$#{ $self->{coords} }) {
    foreach my $atom2 ($atom1+1..$#{ $self->{coords} }) {
      my $string = $self->{elements}->[$atom1] . $self->{elements}->[$atom2];
      if((my $sigma = $rij->{$string}) && (my $epsilon = $eij->{$string})) {
        my $R = $self->distance(atom1 => $atom1, atom2 => $atom2);
        $energy += $epsilon*(($sigma/$R)**12 - ($sigma/$R)**6);
      }
    }
  }
  return $energy;
}


#calculates  LJ-6-12 potential energy based on autodock Rij and Eij parameters
#simply ignores any atom pair involving elements for which parameters are missing (which shouldn't be anything!).
sub LJ_energy_with {
  my ($self, $geo2) = @_;

  my $energy = 0;

  foreach my $atom1 (0..$#{ $self->{coords} }) {
    foreach my $atom2 (0..$#{ $geo2->{coords} }) {
      my $string = $self->{elements}->[$atom1] . $geo2->{elements}->[$atom2];
      if((my $sigma = $rij->{$string}) && (my $epsilon = $eij->{$string})) {
        my $R = $self->distance(atom1 => $atom1, atom2 => $atom2, geometry2 => $geo2);
        $energy += $epsilon*(($sigma/$R)**12 - ($sigma/$R)**6);
      }
    }
  }
  return $energy;
}


#finds minimum energy structure (based on LJ potential) by rotating list of target atoms around bond between $atom1 and $atom2
#TODO: if no list of @targets provided, rotate fragment that starts with atom1!
sub minimize_torsion {
  my $self = shift;
  my %param = @_;
  my ($atom1, $atom2) = ( $param{start_atom}, $param{end_atom} );

  my $targets = $self->get_all_connected($atom1, $atom2);

  my $increment = 5; #angle increment to search over

  my $E_min = 1E20;
  my $angle_min = 0;

  my $axis = $self->get_bond($atom2, $atom1);

  foreach my $count (0..360/$increment + 1) {
    my $angle = $count*$increment;
    $self->center_genrotate($atom1, $axis, deg2rad($increment), $targets);
    my $energy = $self->LJ_energy();
    if($energy < $E_min) {
      $angle_min = $angle;
      $E_min = $energy;
    }
  }
  
  $self->center_genrotate($atom1, $axis, deg2rad($angle_min), $targets);
}


sub get_bond {
    my ($self, $atom1, $atom2, $geometry_2) = @_;
    $geometry_2 //= $self;

    my $pt1 = $self->get_point($atom1);
    my $pt2 = $geometry_2->get_point($atom2);
    
    my $bond = $pt1 - $pt2;
    return $bond;
}


###################
#RMSD related part#
###################
sub RMSD {
    my $self = shift;

    my %params = @_;

    my ($geo2, $heavy_only, 
        $atoms1_ref, $atoms2_ref) = ( $params{ref_geo}, $params{heavy_atoms},
                                      $params{ref_atoms1}, $params{ref_atoms2} );

    $heavy_only //= 0;
    $atoms1_ref //= [0..$#{ $self->{elements} }];
    $atoms2_ref //= [0..$#{ $geo2->{elements} }];

    my $cen1 = $self->get_center($atoms1_ref);
    my $cen2 = $geo2->get_center($atoms2_ref);

    $geo2 = $geo2->copy();

    $self->coord_shift(-1*$cen1);
    $geo2->coord_shift(-1*$cen2);

    for my $i (0..2) {
        map {$atoms1_ref->[$_]->[$i] -= $cen1->[$i]} 
            grep {$atoms1_ref->[$_] !~ /^\d+$/} (0..$#{$atoms1_ref});
        map {$atoms2_ref->[$_]->[$i] -= $cen2->[$i]} 
            grep {$atoms2_ref->[$_] !~ /^\d+$/} (0..$#{$atoms2_ref});
    }

    my $rmsd = $self->_RMSD($geo2, $heavy_only, $atoms1_ref, $atoms2_ref);

    $self->coord_shift($cen2);

    for my $i (0..2) {
        map {$atoms1_ref->[$_]->[$i] += $cen2->[$i]} 
            grep {$atoms1_ref->[$_] !~ /^\d+$/} (0..$#{$atoms1_ref});
    }

    return $rmsd;
}


sub MSD {
    my $self = shift;

    my %params = @_;

    my ($geo2, $heavy_only, 
        $atoms1_ref, $atoms2_ref) = ( $params{ref_geo}, $params{heavy_atoms},
                                      $params{ref_atoms1}, $params{ref_atoms2} );

    my $msd = $self->_RMSD($geo2, $heavy_only, $atoms1_ref, $atoms2_ref, 1);

    return $msd;
}


sub _RMSD {
    my $self = shift;

    my ($geo2, $heavy_only, $atoms1_ref, $atoms2_ref, $no_rot) = @_;

    my $matrix = new Math::MatrixReal(4,4);

    for my $atom (0..$#{$atoms1_ref}) {
        if ($atoms1_ref->[$atom] =~ /^\d+$/ &&
            ($atoms2_ref->[$atom] =~ /^\d+$/) &&
            $heavy_only) {
            if ( ($self->{elements}->[$atoms1_ref->[$atom]] eq 'H')
                && ($geo2->{elements}->[$atoms2_ref->[$atom]] eq 'H') ) {
                next;
            }
        }
        my $pt1 = $atoms1_ref->[$atom] =~ /^\d+$/ ?
                  $self->get_point($atoms1_ref->[$atom]) : $atoms1_ref->[$atom];
        my $pt2 = $atoms2_ref->[$atom] =~ /^\d+$/ ?
                  $geo2->get_point($atoms2_ref->[$atom]) : $atoms2_ref->[$atom];

        $matrix += quat_matrix($pt1, $pt2);
    }

    my ($eigenvalues, $evectors) = $matrix->sym_diagonalize();

    #find smallest of four eigenvalues and save corresponding eigenvectors
    my $sd = 999;
    my $Q = new Math::MatrixReal(1,4);
    foreach my $i (1..4) {
      my $value = $eigenvalues->element($i,1);
      if($value < $sd) {
        $sd = $value;
        $Q = $evectors->column($i);
      }
    }
    
    my $rmsd = 0;
    if($sd > 0) { #to avoid very small negative numbers for sd (-1e-16, etc)
      $rmsd = sqrt($sd/($#{$atoms1_ref}+1));
    }

    my $a = $Q->element(1,1);

    my $w = V($Q->element(2,1), $Q->element(3,1), $Q->element(4,1));

    unless ($no_rot){
        $self->quat_rot($a, $w);
        for my $vec (@$atoms1_ref) {
            if ($vec !~ /^\d+$/) {
                &__point_quat_rot($vec, $a, $w);
            }
        }
    }

    return $rmsd;
}


#This version of RMSD will mirror the molecules with respect to three planes
#Don't use this RMSD unless you know what you are doing, just use normal RMSD.
#This RMSD is particularly designed for non-covalent interaction.
sub RMSD_mirror {
    my $self = shift;

    my %params = @_;

    my ($geo2, $heavy_only, 
        $atoms1_ref, $atoms2_ref) = ( $params{ref_geo}, $params{heavy_atoms},
                                      $params{ref_atoms1}, $params{ref_atoms2} );

    $heavy_only //= 0;
    $atoms1_ref //= [0..$#{ $self->{elements} }];
    $atoms2_ref //= [0..$#{ $geo2->{elements} }];

    my $cen1 = $self->get_center($atoms1_ref);
    my $cen2 = $geo2->get_center($atoms2_ref);

    $self->coord_shift(-1*$cen1);
    $geo2->coord_shift(-1*$cen2);

    #Run regular RMSD comparison as well as mirroring compare_geo across each plane
    #Report minimum value as the RMSD
    my @sd;
    my @Q;

    foreach my $mirrorx (0,1) {
        foreach my $mirrory (0,1) {
            foreach my $mirrorz (0,1) {
                if($mirrorx) {
                    $geo2->mirror_coords('X');
                }
                if($mirrory) {
                    $geo2->mirror_coords('Y');
                }
                if($mirrorz) {
                    $geo2->mirror_coords('Z');
                }
                my $matrix = new Math::MatrixReal(4,4);

                for my $atom (0..$#{$atoms1_ref}) {
                    my $pt1 = $self->get_point($atom);
                    my $pt2 = $geo2->get_point($atom);

                    $matrix += quat_matrix($pt1, $pt2);
                }

                my ($eigenvalues, $evectors) = $matrix->sym_diagonalize();
                #find smallest of four eigenvalues and save corresponding eigenvectors
                my $sd = 999;
                my $Q = new Math::MatrixReal(1,4);
                foreach my $i (1..4) {
                  my $value = $eigenvalues->element($i,1);
                  if($value < $sd) {
                    $sd = $value;
                    $Q = $evectors->column($i);
                  }
                }
                push (@sd, $sd);
                push (@Q, $Q);
            }
        }
    }

    my @idx_sort = sort { $sd[$a] <=> $sd[$b] } 0..$#sd;

    @sd = @sd[@idx_sort];
    @Q = @Q[@idx_sort];
            
    my $rmsd = 0;
    if($sd[0] > 0) { #to avoid very small negative numbers for sd (-1e-16, etc)
       $rmsd = sqrt($sd[0]/($#{$atoms1_ref}+1));
    }else {
       $rmsd = sqrt(-1*$sd[0]/($#{$atoms1_ref}+1));
    }

    my $a = $Q[0]->element(1,1);
    my $w = V($Q[0]->element(2,1), $Q[0]->element(3,1), $Q[0]->element(4,1));

    $self->quat_rot($a, $w);
    
    $self->coord_shift($cen2);

    return $rmsd;
}




#this function map a catalyst to a ts from TS library. The old catalyst
#will be replaced by the new one. Here, $self is the geometry instance for 
#the new catalyst and geo_ref is the TS geometry instance. To replace old
#catalyst, the first atom of the catalyst in the ts should be provided.
#FIXME this is a very priliminary function, and only can be used for some very specific 
#purpose. 
sub map_catalyst {
    my $self = shift;
    my %params = @_;
    my ($geo_ref, $bonds_LJ, $first_cat_atom) = ( $params{geo_ref},
                                                $params{bonds_LJ},
                                                $params{first_cat_atom} );
    
    my $mapped_cata = $self->map_molecule( geo_ref => $params{geo_ref}, 
                                           key_atoms1 => $params{key_atoms1}, 
                                           key_atoms2 => $params{key_atoms2}, 
                                           bonds => $params{bonds} );

    my $num_splice = $#{ $geo_ref->{elements} } - $first_cat_atom + 1;
    my $new_geo = $geo_ref->copy();

    $new_geo->splice_atom($first_cat_atom, $num_splice);
    $new_geo->append($mapped_cata->copy());
    #FIXME the catatlyst rotatation was removed 
    for my $bond_LJ (@$bonds_LJ) {
        my @bond_LJ = map {$_ + $first_cat_atom} @$bond_LJ;
        $new_geo->minimize_torsion(@bond_LJ);
    }
    return $new_geo;
}


sub get_center {
    my ($self, $groups) = @_;

    $groups //= [0.. $#{ $self->{elements} }];

    my @xyz_groups = grep { $_ !~ /^\d+$/} @$groups;
    my @groups = grep{ $_ =~ /^\d+$/ } @$groups;

    my $COM = V(0, 0, 0);
    for my $atom (@groups) {$COM += $self->{coords}->[$atom];}

    for my $point (@xyz_groups) {$COM += $point;}

    $COM /= $#{ $groups } + 1;

    return $COM;
}


#Find shortest path between atom1 and atom2
#Performs a breadth first search, returns length of shortest path
#returns -1 if no path found
sub shortest_path {
  my ($self, $atom1, $atom2) = @_;
  my @positions = ($atom1);
  my %visited = ( $atom1 => '0') ;	#keys are numbers of the visited atoms, values are the corresponding depths

  #loop until @positions is empty
  while(@positions) {
    my $position = shift(@positions);
    if ($position eq $atom2) {
      return $visited{$position};
    }
    foreach my $atom (@{$self->{connection}->[$position]}) { 	#if not found yet, grab all atoms connected to current atom and add to queue (unless already visited)
      if(!exists $visited{$atom}) {	#skip over element, just add atom numbers
        push(@positions, $atom);
        $visited{$atom} = $visited{$position} + 1;
      }
    }
  }
  return -1;	#return -1 if no path found in BFS
} #end shortest_path


sub rotatable_bonds {
    my $self = shift;

    my ($ref1, $ref2) = @_;

    my %bonds;
    my $empty = 1;

    for my $activei (@$ref1) {
        for my $activej (@$ref2) {
            my @path = ({$activei => -1});
            while (@path) {
                my @newpath = ();
                for my $path (@path) {
                    my ($head) = grep {$path->{$_} < 0} keys %{ $path };
                    for my $atom_next (@{$self->{connection}->[$head]}) {
                        my $new_path = { %$path };
                        $new_path->{$head} = $atom_next;

                        if ($atom_next == $activej) {
                            if ($empty) {
                                my @keys = keys %{ $new_path };
                                @bonds{@keys} = map {$new_path->{$_}} @keys;
                                $empty = 0;
                            }else {
                                for my $key (keys %bonds) {
                                    if (exists $new_path->{$key} &&
                                        ($new_path->{$key} == $bonds{$key})) {
                                            next;
                                    }else {
                                        delete $bonds{$key};
                                    }
                                }
                            }
                        }elsif (! exists $new_path->{$atom_next}) { 
                            $new_path->{$atom_next} = -1;
                            push (@newpath, $new_path);
                        }
                    }
                }
                @path = @newpath;
            }
        }
    }

    for my $key (keys %bonds) {
        my $ele1 = $self->{elements}->[$key];
        my $ele2 = $self->{elements}->[$bonds{$key}];
        unless (@{$self->{connection}->[$key]} >= $CONNECTIVITY->{$ele1} ||
                @{$self->{connection}->[$bonds{$key}]} >= $CONNECTIVITY->{$ele2}) {
            my $d = $self->distance(atom1 => $key, atom2 => $bonds{$key});
            my $d_ref = $UNROTATABLE_BOND->{$ele1.$ele2} ?
                        $UNROTATABLE_BOND->{$ele1.$ele2} : $UNROTATABLE_BOND->{$ele2.$ele1} ?
                                                           $UNROTATABLE_BOND->{$ele2.$ele1} : 0;
            if ($d_ref) {
                if ($d < $d_ref) {
                    delete $bonds{$key};
                }
            }else {
                my $msg = "Unrotatable bond criterion $ele1-$ele2 is not implemented, " .
                          "This bond will be viewed as rotatable bond\n";
                warn($msg);
            }
        }
    }
    return \%bonds;
}


sub bare_backbone {
    my $self = shift;

    my ($active_centers) = @_;

    my %backbone;
    my @active_centers = @$active_centers;

    while (@active_centers) {
        my $activei = shift @active_centers;
        my $activej = $active_centers[0] || $activei;

        my @path = ({$activei => 1, head => $activei});
        while (@path) {
            my @newpath = ();
            for my $path (@path) {
                for my $atom_next (@{$self->{connection}->[$path->{head}]}) {
                    my $new_path = { %$path };
                    $new_path->{$atom_next} = 1;

                    if (($atom_next == $activei || ($atom_next == $activej)) &&
                        (keys %{ $new_path } > 3)) {
                        @backbone{keys %{ $new_path }} = ();
                    }elsif (! exists $path->{$atom_next}) { 

                        $new_path->{head} = $atom_next;
                        push (@newpath, $new_path);
                    }
                }
            }
            @path = @newpath;
        }
    }

    delete $backbone{head};

    %backbone = map { $_ => 1 } @$active_centers if (keys %backbone == 0);

    return \%backbone;
}


sub printXYZ {
    my ($self, $filename, $comment) = @_;
    $comment //= $self->{name};
    my $num_atoms = $#{ $self->{elements} } + 1;

    my $handle;
    if($filename) {
        open $handle, ">$filename" or die "Can't open $filename\n";
    }else {
        $handle = *STDOUT;
    }
    print $handle "$num_atoms\n$comment\n";
    foreach my $atom (0..$#{ $self->{elements} }) {
        printf $handle "%2s%14.6f%14.6f%14.6f\n", ($self->{elements}->[$atom], @{ $self->{coords}->[$atom] });
    }
    close $handle if $filename;
}


#Writes com file
#write_com(route, comment, charge, multiplicity, ref_to_coords, footer, flag)
#write_com(route, comment, charge, multiplicity, ref_to_coords, footer, flag, filename)
#where footer contains anything that goes after the coords (gen basis specification, mod redundant commands, etc)
#flag = 0 will print only elements and coords
#flag = 1 will print 0/-1 column as well
sub write_com {
    my $self = shift;
    my %params = @_;
    my ($comment, $route, $charge, $mult, $footer, $flag, $filename) = ( $params{comment},
                                                                         $params{route},
                                                                         $params{charge},
                                                                         $params{mult},
                                                                         $params{footer},
                                                                         $params{print_flag},
                                                                         $params{filename} );
    my $fh;
    $filename && open ($fh, ">$filename") || ($fh = *STDOUT);
    
    print $fh "$route\n\n";
    print $fh "$comment\n\n";
    print $fh "$charge $mult\n";

    foreach my $atom (0..$#{ $self->{elements} }) {
        if ($flag) {
            printf $fh "%-2s%4s%14.6f%14.6f%14.6f\n", ($self->{elements}->[$atom], 
                                                       $self->{flags}->[$atom],
                                                       @{ $self->{coords}->[$atom] });
        }else { 
            printf $fh "%2s%14.6f%14.6f%14.6f\n", ($self->{elements}->[$atom], 
                                                       @{ $self->{coords}->[$atom] });
        }
    }

    print $fh "\n";
    if($footer) {
      print $fh "$footer\n";
    }
    print $fh "\n\n";

    close ($fh);
}


sub flatten {
    my ($self) = @_;
    
    my $num_atoms = $#{ $self->{elements} } + 1;
    my $geometry = "$num_atoms\\\\n\\\\n";
    foreach my $atom (0..$#{ $self->{elements} }) {
        $geometry .= "$self->{elements}->[$atom] $self->{coords}->[$atom]->[0] ".
                     "$self->{coords}->[$atom]->[1] $self->{coords}->[$atom]->[2]\\\\n";
    }

    return $geometry;
}




sub union {
    my ($a, $b) = @_;
    $a //= [];
    $b //= [];

    my @union = ();
    my %union = ();

    for my $e (@$a) { $union{$e} = 1 }
    for my $e (@$b) { $union{$e} = 1 }

    return [sort keys %union];
}


sub quat_matrix {
    my ($pt1, $pt2) = @_;
    my ($xm, $ym, $zm) = @{ $pt1 - $pt2 };
    my ($xp, $yp, $zp) = @{ $pt1 + $pt2 };

    my $temp_matrix = Math::MatrixReal->new_from_rows( 
        [[$xm*$xm + $ym*$ym + $zm*$zm, $yp*$zm - $ym*$zp,          $xm*$zp - $xp*$zm,           $xp*$ym - $xm*$yp],
        [$yp*$zm - $ym*$zp,           $yp*$yp + $zp*$zp + $xm*$xm,$xm*$ym - $xp*$yp,           $xm*$zm - $xp*$zp],
        [$xm*$zp - $xp*$zm,           $xm*$ym - $xp*$yp,          $xp*$xp + $zp*$zp + $ym*$ym, $ym*$zm - $yp*$zp], 
        [$xp*$ym - $xm*$yp,           $xm*$zm - $xp*$zp,          $ym*$zm - $yp*$zp,           $xp*$xp + $yp*$yp + $zm*$zm]]
    );

    return $temp_matrix;
}


#########################################
#some internal function you should never#
#call outside                           #
#########################################

sub __point_quat_rot {
    my ($vec, $a, $w) = @_;

    my $wx = $w x $vec;
    my $new_vec = $vec + 2*$a*$wx + 2*($w x $wx);
    #This is to maintain the array ref
    map {$vec->[$_] = $new_vec->[$_]} (0..2);
}    


package AaronTools::NanoTube;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
our @ISA = qw(AaronTools::Geometry);

my $CC = 1.415;
my $CH = 1.08;

sub new {
    #initiate
    my $class = shift;
    my %params = @_;
    my $name = $params{name} ? 
               "$params{name}-$params{width}-$params{length}" : 
               "nt-$params{width}-$params{length}";

    my $self = new AaronTools::Geometry(name => $name);
   
    $self->{width} = $params{width};
    $self->{length} = $params{length};
    $self->{radius} = $params{radius};
    $self->{angular_offset} = $params{angular_offset} // 0; 

    my $fragment = $self->{radius} ? 1 : 0;
    $self->{radius} //= ($self->{width} >= 2) ? 
                        newton($CC, $self->{width}) : 0;
    unless ($self->{radius}) {die ("Can't build nanotube smaller than (2,2)!\n");}
        
    bless $self, $class;

    #make new atom geometry 
    my $atom = new AaronTools::Geometry(name => 'carbon',
                                        elements => ['C'],
                                        coords => [[$self->{radius}, 0, 0]]);

    #build nano tube;
    my $CC_angle = 2*asin($CC/(2*$self->radius()));
    my $CC_halfangle = 2*asin($CC/(4*$self->radius()));
    my $CC_side = $CC*sqrt(3.0)/2.0;

    my $a = $self->angular_offset();
    my $angle = -($self->width()/2+$self->width()-1)*$CC_angle - 
                $self->angular_offset()*2*($CC_angle+$CC_halfangle);
    $atom->rotate('z', $angle);

    my $shift = V(0, 0, $CC_side*($self->length()-1)/2);
    $atom->coord_shift($shift);

    my $angle_tally = 0;


    for(my $row=0; $row<$self->length(); $row++) {
        if($row%2==1) {
            if($row!=$self->length()-1 || $fragment==0) {
                $self->append($atom->copy());
            }
            $atom->rotate('z', $CC_angle+2*$CC_halfangle);
            $angle_tally += $CC_angle+2*$CC_halfangle;
        } else {
            $atom->rotate('z', $CC_halfangle);
            $angle_tally += $CC_halfangle;
        }
        for (my $ring=0; $ring<$self->width(); $ring++) {
            if($row!=$self->length()-1 || $row%2!=1 || $ring != $self->width()-1 || $fragment==0) {
                $self->append($atom->copy());
            }
            $atom->rotate('z', $CC_angle);
            $angle_tally += $CC_angle;
            if($row%2!=1 || $ring != $self->width()-1) {
                $self->append($atom->copy());
            }
            $atom->rotate('z', $CC_angle+2*$CC_halfangle);
            $angle_tally += $CC_angle+2*$CC_halfangle;
        }
    
        #Reset and shift
    #    slide(-$CC_side);
        $atom->coord_shift(V(0, 0, -$CC_side)); 
        $atom->rotate('z', -$angle_tally);
        $angle_tally = 0;
    }

    #Cap open valences
    my $Hatom = new AaronTools::Geometry(name => 'hydrogen',
                                        elements => ['H'],
                                        flags => [0],
                                        coords => [[0, 0, 0]]);
    
    my $numCs = $#{ $self->{elements} };
    my $Hatoms = new AaronTools::Geometry(name => 'Hs');
    foreach my $atom1 (0..$numCs) {
        my $vector = V(0, 0, 0);
        my $neighbors = 0;
        foreach my $atom2 (0..$numCs) {
            if($atom1 != $atom2) {
                if($self->distance(atom1 => $atom1, atom2 => $atom2) < $CC+0.1) {
                    $neighbors++;
                    $vector += $self->get_bond($atom2, $atom1);
                }
            }
        }
        if($neighbors < 3) {
            my $norm = abs($vector);
            $vector /= $norm;
            my $coord = [ @{ $self->get_point($atom1) - $vector } ];
            $Hatom->update_coords(targets => [0], coords => [$coord]); 
            $Hatoms->append($Hatom->copy());
        }
        if($neighbors < 2) {
            die "Dangling carbon $atom1!";
        }
        if($neighbors > 4) {
            die "Too many neighbors for atom $atom1 (radius too small to accommodate $self->width() rings)";
        }
    }
    $self->append($Hatoms->copy());
    $self->refresh_connected();
    return $self;
}


sub copy {
    my $self = shift;
    my $new =  new AaronTools::Geometry( name => $self->{name},
                                         elements => [ @{ $self->{elements} } ],
                                         flags => [ @{ $self->{flags} }],
                                         coords => [ map { [ @$_ ] } @{ $self->{coords} } ],
                                         connection => [ map { [ @$_ ] } @{ $self->{connection} } ],
                                         constraints => [ map { [ @$_ ] } @{ $self->{constraints} } ] );
    for my $key ("width", "length", "radius", "angular_offset") {
        $new->{$key} = $self->{$key};
    }
    bless $new, "AaronTools::NanoTube";
    return $new;
};


sub width {
    my $self = shift;
    return $self->{width};
}


sub length {
    my $self = shift;
    return $self->{length};
}


sub radius {
    my $self = shift;
    return $self->{radius};
}


sub angular_offset {
    my $self = shift;
    return $self->{angular_offset};
}


#Dirty Newton solver to get radius of closed CNTs
sub newton {
  my ($CC, $width) = @_;
  #Threshold for Newton-Raphson solver to get radius
  my $THRESH = 1E-10;

  my $lastradius = 3*$CC*$width/(2*pi()); #starting guess from conventional formula
  my $old_gap = get_CNT_gap($lastradius, $CC, $width);
  my $radius = $lastradius + 0.01; #arbitrary step to get process started
  my $gap = get_CNT_gap($radius, $CC, $width);

  #Simple Newton solver using very crude finite difference derivative
  while(abs($gap) > $THRESH) {
    my $newradius = $radius - $gap*($radius - $lastradius)/($gap - $old_gap);
    $old_gap = $gap;
    $gap = get_CNT_gap($newradius, $CC, $width);
    $lastradius = $radius;
    $radius = $newradius;
  }
  return $radius;
}


#Function to be minimized to get the radius of closed CNT
sub get_CNT_gap {
  my ($guess, $CC, $width) = @_;
  my $value = asin($CC/(2*$guess)) + asin($CC/(4*$guess)) - pi()/(2*$width);
  return $value;
}



package AaronTools::Substituent;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
our @ISA = qw(AaronTools::Geometry);

sub new {
    my $class = shift;
    my %params = @_;

    my $self = new AaronTools::Geometry();
    delete $self->{constraints};
    bless $self, $class;

    if (exists $params{name}) {
        $self->set_name($params{name});
        if (-f "$AARON/Subs/$self->{name}.xyz") {
            $self->read_geometry("$AARON/Subs/$self->{name}.xyz");
        }
    }

    $self->{end} = $params{end};
    #This is the potential substituent in the furture
    $self->{sub} = $params{sub};

    return $self;
}


sub copy {
    my $self = shift;

    my $new =  new AaronTools::Geometry( name => $self->{name},
                                         elements => [ @{ $self->{elements} } ],
                                         flags => [ @{ $self->{flags} }],
                                         coords => [ map { [ @$_ ] } @{ $self->{coords} } ],
                                         connection => [ map { [ @$_ ] } @{ $self->{connection} } ] );

    bless $new, "AaronTools::Substituent";
    $new->{end} = $self->{end};
    $new->{sub} = $self->{sub};
    $new->{conformer_num} = $self->{conformer_num};
    $new->{conformer_angle} = $self->{conformer_angle};
    return $new;
};


sub end {
    my $self = shift;
    return $self->{end};
}


sub read_geometry {
    my ($self, $file) = @_;

    my ($elements, $flags, $coords, $constraints, $ligand,
        $TM, $key_atoms, $bonds, $conformer) = AaronTools::FileReader::grab_coords($file);

    $self->{elements} = $elements;
    $self->{flags} = $flags;
    $self->{coords} = $coords;
    $self->{conformer_num} = $conformer->[0];
    $self->{conformer_angle} = $conformer->[1];

    $self->refresh_connected();
}


sub compare_lib {
    my $self = shift;

    my $subs = {};

    open (my $fh, "<$AARON/Subs/subs") or die "Cannot open AARON/Subs/subs";

    while (<$fh>) {
        chomp;
        if ($_ =~ /[0-9a-zA-Z]/) {
            my $name = $_;
            my $sub = new AaronTools::Substituent( name => $name );
            delete $sub->{coords};
            delete $sub->{flags};
            $subs->{$name} = $sub;
        }
    }

    for my $sub (keys %{ $subs }) {
        if ($#{ $subs->{$sub}->{elements} } != $#{ $self->{elements} }) {
            delete $subs->{$sub};
        }else {
            $subs->{$sub}->{visited} = {};
            $subs->{$sub}->{open_set} = { 0 => $subs->{$sub}->{elements}->[0] };
        }
    }
    #initiate $self_sub
    my $self_sub = {};
    $self_sub->{connection} = $self->{connection};
    $self_sub->{elements} = $self->{elements};
    $self_sub->{visited} = {};
    $self_sub->{open_set} = { 0 => $self->{elements}->[0] };

    while (%{ $subs } && %{ $self_sub->{open_set} }) {
        for my $sub (keys %{ $subs }) {
            unless(&_same_nodes($subs->{$sub}->{open_set}, $self_sub->{open_set})) {
                delete $subs->{$sub};
                next;
            }else {
                &_move_to_next_layer($subs->{$sub});
            }
        }
        &_move_to_next_layer($self_sub);
    }

    my @subs = keys %{ $subs };
    if ($#subs < 0) {
        print "Cannot detect what kind of sunstituent it is, no conformer information got.\n";
    }elsif ($#subs == 0) {
       $self->{name} = $subs[0];
       $self->{conformer_num} = $subs->{$subs[0]}->{conformer_num};
       $self->{conformer_angle} = $subs->{$subs[0]}->{conformer_angle};
    }else {
       print "Multiple similar substituent was found, but can not tell which one it is (No conformer information): \n";
       print "@subs\n";
    }
}


#align the substituent to a bond of geometry object
sub _align_on_geometry {
    my $self = shift;
    my %params = @_;

    my ($geo, $target, $end) = ($params{geo}, $params{target}, $params{end});

    #Rotate to align along nearest_neighbor-target bond then shift sub_coords to nearest_neighbor position
    my $nearst_v = $geo->get_point($end);
    my $target_v = $geo->get_point($target);

    my $bond_axis = $target_v - $nearst_v;

    $bond_axis /= $bond_axis->norm();

    #sub_coords are aligned along x-axis, so find rotation axis that transforms x-axis to bond_axis
    my $v_x = V(1,0,0);
    my $cross = $v_x x $bond_axis;
    my $angle = atan2($bond_axis, $v_x);
    
    $self->genrotate($cross, $angle);
    $self->coord_shift($nearst_v);

    my $current_distance = $self->distance( atom1 => 0, atom2 => $end,
                                           geometry2 => $geo);
    my $current_bond = $self->get_bond( 0, $end, $geo);

    my $new_distance = $radii->{$self->{elements}->[0]} 
                       + $radii->{$geo->{elements}->[$end]};

    my $difference = $new_distance - $current_distance;

    my $v12 = $current_bond * $difference / $new_distance;

    $self->coord_shift($v12);
}


#######################
#Some useful functions#
#######################
sub _same_nodes {
    my ($set1, $set2) = @_;

    my %set1 = %{ $set1 };
    my %set2 = %{ $set2 };

    my $same_nodes = 1;

    if (keys %set1 != keys %set2) {
        $same_nodes = 0;
    }else {
        my @keys_set1 = sort {$a <=> $b} keys %set1;
        while(%set1 && $same_nodes) {
            my $i = shift @keys_set1;
            my $j;
            my $found = 0;
            for my $j_temp (keys %set2) {
                if ($set1{$i} eq $set2{$j_temp}) {$found = 1; $j=$j_temp; last;}
            }
            unless ($found) {
                $same_nodes = 0
            }else{ delete $set1{$i}; delete $set2{$j}; }
        }
    }

    if (%set2) {$same_nodes = 0};

    return $same_nodes;
}


sub _move_to_next_layer {
    my ($set) = @_;

    for my $key (keys %{ $set->{open_set} } ) {
        $set->{visited}->{$key} = $set->{open_set}->{$key};
        delete $set->{open_set}->{$key};

        for my $connected (@{ $set->{connection}->[$key] }) {
            if (! exists $set->{visited}->{$connected}) {
                $set->{open_set}->{$connected} = $set->{elements}->[$connected];
            }
        }
    }
}
