#Contributors: Yanfei Guan and Steven Wheeler

use lib $ENV{'QCHASM'};
use lib $ENV{'PERL_LIB'};

use AaronTools::Constants qw(CUTOFF NAMES);
use AaronTools::Atoms qw (:BASIC :LJ);

my $QCHASM = $ENV{'QCHASM'};
my $mass = MASS;
my $CUTOFF = CUTOFF;
my $TMETAL = TMETAL;
my $RADII = RADII;
my $CONNECTIVITY = CONNECTIVITY;
my $radii = RADII;

my $SUBSTRATE = NAMES->{SUBSTRATE};
my $LIGAND = NAMES->{LIGAND};
my $CENTER = NAMES->{CENTER};


package AaronTools::Catalysis;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
use Math::MatrixReal;
use AaronTools::Geometry;
use Data::Dumper;
our @ISA = qw(AaronTools::Geometry);

sub new {
    my $class =shift;
    my %params = @_;

    my $substituents = $params{substituents};

    my $self = new AaronTools::Geometry();
    bless $self, $class;

    if (exists $params{name}) {
        $self->set_name($params{name});
        if (-f "$self->{name}.xyz") {
            $self->read_geometry("$self->{name}.xyz");
            #ligand backbone and sub
            $self->detect_component();
            $self->ligand()->set_substituents($substituents->{ligand});
            $self->ligand()->detect_backbone_subs( no_new_subs => $params{no_new_subs} );
            #substrate subs

            $self->substrate()->set_substituents($substituents->{substrate});
            $self->substrate()->detect_backbone_subs();
        }
    }

    $self->{conf_num} = 1;

    $self->rebuild();

    $self->refresh_connected();

    return $self;
}


#This is to update the catalysis from a file like .xyz or .log
#This works better with .xyz file with comment specifing the information
#of the system. For .log file or .com file, if the atom order is with ligand
#first then TM then substrate, it will work fine. But if the order is not like
#this,the catalysis will be no different from a general geometry object
#F:12-19 L:0-23 T:24 K:0,3
#This is especially used to update ligand, metal(center) and substrate from catalysis
#coordinates
sub detect_component {
    my $self = shift;

    my $TM;
    for my $atomi (0..$#{ $self->{elements} }) {
        if (exists $TMETAL->{$self->{elements}->[$atomi]}) {
            $TM = $atomi;
        }
    }

    if ($self->{ligand_atoms}) {

        if (! $self->{center_atom}) {
            if ($TM) {
                $self->{center_atom} = $TM;
                $self->{substrate_atoms} = [0..$TM-1];
            }else {
                print "No transition metal was found in the geometry. ".
                      "Catalysis is found to be a pure organic system. " .
                      "If this is a Si, P or other non-metal atom centered ".
                      "system, please specify that in the .xyz file.\n";
                $self->{substrate_atoms} = [0..$self->{ligand_atoms}->[0]-1];
            }
        }else {
            $self->{substrate_atoms} = [0..$self->{center_atom}-1];
        }
    }else {
        if (! $self->{center_atom}) {
            if ($TM) {
                $self->{center_atom} = $TM;
                $self->{ligand_atoms} = [$TM+1..$#{ $self->{elements} }];
                $self->{substrate_atoms} = [0..$TM-1];
            }else {
                die "This is a pure organic system " .
                    "You must specify the ligand atoms in the xyz file";

            }
        }else {
            $self->{ligand_atoms} = [$self->{center_atom}+1..$#{ $self->{elements} }];
            $self->{substrate_atoms} = [0..$self->{center_atom}-1];
        }
    }

    #make center attribute
    if ($self->{center_atom}) {
        $self->{center} = $self->subgeo([$self->{center_atom}]);
        bless $self->{center}, "AaronTools::Component";
    }

    #make ligand attribute
    $self->{ligand} = $self->subgeo($self->{ligand_atoms});
    bless $self->{ligand}, "AaronTools::Ligand";
    $self->{ligand}->refresh_connected();

    if ($self->{ligand_keyatoms}) {
        $self->{ligand}->{active_centers} = $self->{ligand_keyatoms};
    }elsif ($self->{center_atom}) {
        my @key_atoms;

        for my $i (0..$#{ $self->{ligand}->{elements} }) {

            my $d = $self->ligand()->distance(atom1 => $i,
                                              atom2 => 0,
                                              geometry2 => $self->{center});

            if ($d < $RADII->{$self->{ligand}->{elements}->[$i]} +
                     $RADII->{$self->{center}->{elements}->[0]} +
                     $CUTOFF->{D_CUTOFF}) {
                push (@key_atoms, $i);
            }
        }

        unless (@key_atoms) {
            die "No connection detected between ligand and metal " .
                "You probably need to write the key atoms in the xyz file.\n";
        }else {
            $self->{ligand_keyatoms} = [[@key_atoms]];
            $self->{ligand}->{active_centers} = $self->{ligand_keyatoms};
        }
    }else {
        die "For a poure organic system, you are required to provide " .
            "active center information for the catalyst";
    }


    #make substrate attribute
    $self->{substrate} = $self->subgeo($self->{substrate_atoms});
    bless $self->{substrate}, "AaronTools::Substrate";

    $self->{substrate}->refresh_connected();

    #put constraint in each component
    if ($self->{constraints}) {
        my @ligand_con;
        my @substrate_con;

        for my $bond_d (@{ $self->{constraints} }) {
            #all substrate and ligand exist;
            my $bond = $bond_d->[0];
            my $d = $bond_d->[1];

            if (grep (/^$bond->[0]$/, @{ $self->{substrate_atoms} }) &&
                grep (/^$bond->[1]$/, @{ $self->{substrate_atoms} })) {
                push (@substrate_con, [$bond, $d]);
            }elsif (grep (/^$bond->[0]$/, @{ $self->{ligand_atoms} }) &&
                    grep (/^$bond->[1]$/, @{ $self->{ligand_atoms} })) {
                    push (@ligand_con,
                          [[ map { $_ - $self->{ligand_atoms}->[0] } @$bond ], $d]);
            }elsif (grep (/^$bond->[0]$/, @{ $self->{ligand_atoms} }) &&
                    grep (/^$bond->[1]$/, @{ $self->{substrate_atoms} })) {
                    push (@ligand_con,
                          [[$bond->[0] - $self->{ligand_atoms}->[0], 'S'], $d]);
                    push (@substrate_con, [[$bond->[1], 'L'], $d]);
            }elsif (grep (/^$bond->[0]$/, @{ $self->{substrate_atoms} }) &&
                    grep (/^$bond->[1]$/, @{ $self->{ligand_atoms} })) {
                    push (@ligand_con,
                          [[$bond->[1] - $self->{ligand_atoms}->[0], 'S'], $d]);
                    push (@substrate_con, [[$bond->[0], 'L'], $d]);
            }elsif (grep (/^$bond->[0]$/, @{ $self->{substrate_atoms} }) &&
                    grep (/^$bond->[1]$/, ($self->{center_atom}))) {
                    push (@substrate_con, [[$bond->[0], 'C'], $d]);
            }elsif (grep (/^$bond->[1]$/, @{ $self->{substrate_atoms} }) &&
                    grep (/^$bond->[0]$/, ($self->{center_atom}))) {
                    push (@substrate_con, [[$bond->[1], 'C'], $d]);
            }elsif (grep (/^$bond->[0]$/, @{ $self->{ligand_atoms} }) &&
                    grep (/^$bond->[1]$/, ($self->{center_atom}))) {
                    push (@ligand_con,
                          [[$bond->[0] - $self->{ligand_atoms}->[0], 'C'], $d]);
            }elsif (grep (/^$bond->[1]$/, @{ $self->{ligand_atoms} }) &&
                    grep (/^$bond->[0]$/, ($self->{center_atom}))) {
                    push (@ligand_con,
                          [[$bond->[1] - $self->{ligand_atoms}->[0], 'C'], $d]);
            }
        }

        if ($self->ligand()) {
            $self->ligand()->{constraints} = [@ligand_con];
            $self->ligand()->{substituents} = {};
        }

        if ($self->substrate()) {
            $self->substrate()->{constraints} = [@substrate_con];
            $self->substrate()->{substituents} = {};
        }

    }
}


#This function is opposite to the@split above, it is used to combine the

#ligand, metal and substrate
sub rebuild {
    my $self = shift;


    $self->rebuild_coords();

    $self->refresh_connected();

    $self->{substrate_atoms} = [0..$#{ $self->{substrate}->{elements} }];
    if ($self->center()) {
        $self->{center_atom} = $#{ $self->{substrate}->{elements} } + 1;
        $self->{ligand_atoms} = [$self->{center_atom}+1..$#{ $self->{elements} }];
    }else {
        $self->{ligand_atoms} = [$#{ $self->{substrate}->{elements} } + 1..
                                 $#{ $self->{elements} }];
    }

    my @ligand_con = @{ $self->ligand()->{constraints} }
                        if $self->ligand()->{constraints};
    my @substrate_con = @{ $self->substrate()->{constraints} }
                            if $self->substrate()->{constraints};

    my @ligand_con_sub = grep { $_->[0]->[1] eq 'S' } @ligand_con;
    my @substrate_con_li = grep { $_->[0]->[1] eq 'L' } @substrate_con;
    my @ligand_con_center = grep { $_->[0]->[1] eq 'C' } @ligand_con;
    my @substrate_con_center = grep { $_->[0]->[1] eq 'C' } @substrate_con;

    my @constraint;

    if ($#ligand_con_sub == $#substrate_con_li) {
        my @temp = map { [[$ligand_con_sub[$_][0][0] + $self->{ligand_atoms}->[0],
                          $substrate_con_li[$_][0][0]], $ligand_con_sub[$_][1]] }
                        (0..$#ligand_con_sub);
        push (@constraint, @temp);
    }

    push (@constraint,
          map {[[$_->[0]->[0] + $self->{ligand_atoms}->[0], $self->{center_atom}], $_->[1]]}
            @ligand_con_center);

    push (@constraint,
          map {[[$_->[0]->[0], $self->{center_atom}], $_->[1]]} @substrate_con_center);

    my @li_li_con = grep { $_->[0]->[1] =~ /^\d+$/ } @ligand_con;
    my @sub_sub_con = grep { $_->[0]->[1] =~ /^\d+$/ } @substrate_con;

    push (@constraint, map { [[$_->[0]->[0] + $self->{ligand_atoms}->[0],
                              $_->[0]->[1] + $self->{ligand_atoms}->[0]], $_->[1]] }
                           @li_li_con);

    push (@constraint, @sub_sub_con);

    $self->{constraints} = [@constraint];
}


sub rebuild_coords {

    my $self = shift;
    $self->ligand()->rebuild();
    $self->substrate()->rebuild();

    my @elements = @{ $self->substrate()->{elements} };
    my @flags = @{ $self->substrate()->{flags} };
    my @coords = @{ $self->substrate()->{coords} };

    if ($self->center()) {
        push (@elements, @{ $self->center()->{elements} });
        push (@flags, @{ $self->center()->{flags} });
        push (@coords, @{ $self->center()->{coords} });
    }

    push (@elements, @{ $self->ligand()->{elements} });
    push (@flags, @{ $self->ligand()->{flags} });
    push (@coords, @{ $self->ligand()->{coords} });

    $self->{elements} = [ @elements ];
    $self->{flags} = [ @flags ];
    $self->{coords} = [ @coords ];

    $self->refresh_connected();
}


sub copy {
    my $self = shift;

    my $new =  new AaronTools::Geometry();
    bless $new, "AaronTools::Catalysis";

    $new->{name} = $self->{name};
    $new->{constraints} = [ map { [ @$_ ] } @{ $self->{constraints} } ];

    $new->{ligand} = $self->{ligand}->copy();
    $new->{center} = $self->{center}->copy() if $self->{center};
    $new->{substrate} = $self->{substrate}->copy();

    $new->{ligand_atoms} = [@{ $self->{ligand_atoms} }];
    $new->{center_atom} = $self->{center_atom};
    $new->{substrate_atoms} = [@{ $self->{substrate_atoms} }];

    $new->{ligand_keyatoms} = [map { [@$_] } @{ $self->{ligand_keyatoms} }];
    $new->rebuild();
    $new->refresh_connected();

    return $new;
}


sub type {
    my $self = shift;
    if ($self->{center}) {
        return 'TM';
    }else {
        return 'Org';
    }
}


sub ligand {
    my $self = shift;
    return $self->{ligand};
}


sub center {
    my $self = shift;
    return $self->{center};
}


sub substrate {
    my $self = shift;
    return $self->{substrate};
}


sub read_geometry {
    my $self = shift;
    my ($file) = @_;
    my ($elements, $flags, $coords, $constraints, $ligand,
        $TM, $key_atoms) = AaronTools::FileReader::grab_coords($file);

    $self->{elements} = $elements if $elements;
    $self->{flags} = $flags if $flags;
    $self->{coords} = $coords if $coords;
    if ($constraints) {
        my @constraints;
        for my $constraint (@$constraints) {
            my $d = $self->distance( atom1 => $constraint->[0],
                                     atom2 => $constraint->[1] );
            push(@constraints, [$constraint, $d]);
            $self->{constraints} = [@constraints];
        }
    }
    $self->{ligand_atoms} = $ligand if $ligand;
    $self->{center_atom} = $TM if $TM;
    $self->{ligand_keyatoms} = $key_atoms if $key_atoms;

    $self->refresh_connected();
}


#This is to update the geometry to each componenet as well as the backbone
#and substituents
sub update_geometry {
    my $self = shift;
    my ($file) = @_;

    $self->read_geometry($file);

    $self->_update_geometry();

}


sub conformer_geometry {
    my $self = shift;
    my ($conf) = @_;

    $self->{coords} = $conf->{coords};
    $self->_update_geometry();
}


sub _update_geometry {
    my $self = shift;

    $self->substrate()->{coords} =
        [ @{$self->{coords}}[@{ $self->{substrate_atoms} }] ];

    $self->center()->{coords} =
        [ ${$self->{coords}}[$self->{center_atom}] ];

    $self->ligand()->{coords} =
        [ @{$self->{coords}}[@{ $self->{ligand_atoms} }] ];

    $self->substrate()->update_geometry();
    $self->ligand()->update_geometry();
}


#set all flags to be 0
sub freeze_all_atoms {
    my $self = shift;

    $self->substrate()->freeze_all_atoms();
    $self->ligand()->freeze_all_atoms();
    $self->center()->{flags} = [-1];

    $self->{flags} = [ map { -1 } @{ $self->{flags} } ];
}


#this is to replace the ligand of the catalysis with a new ligand
#The new ligand will be map to the old one by the key atoms.
#If there are two key atoms, after mapping, the ligand will be rotated
#to minize the rmsd to find the best position
#NOTE: only use for centered catalysis system
#FIXME implement case with three or more key atoms.
sub map_ligand {
    my $self = shift;

    my ($ligand) = @_;

    $self->_map_ligand($ligand);

    $self->remove_clash();

    $self->rebuild()
    #FIXME more atoms case implemented in the furture
}


#FIXME this is only a priliminary version. You should only use this function
#in the case you know what you are doing.
sub _map_ligand {
    my $self = shift;

    my ($ligand_ref) = $self->ligand()->copy();
    my ($ligand) = @_;

    my @key_atoms1 = map { [@$_] } @{ $ligand->{active_centers} };
    my @key_atoms2 = map { [@$_] } @{ $ligand_ref->{active_centers} };

    my @key_total1 = map { @$_ } @key_atoms1;
    my @key_total2 = map { @$_ } @key_atoms2;

    my $i=0;

    unless (@key_atoms1 == @key_atoms2) {
        my $msg = "number of active center blocks for ligand mapping " .
                  "doesn't match." .
                  "The active center will be chunked.\n";
        if (@key_atoms1 < @key_atoms2) {
            @key_atoms2 = @key_atoms2[0..$#key_atoms1];
        }else {
            @key_atoms1 = @key_atoms1[0..$#key_atoms2];
        }
    }

    my $ref1_center = shift @key_atoms1;
    my $ref2_center = shift @key_atoms2;
    if (@$ref1_center != @$ref2_center) {
        warn("reference atom numer not equal for ligand map!\n");
        next;
    }

    $ref1_center = $self->_map_center($ligand, $ref1_center, $ref2_center);

    if (@key_atoms1) {
        $ligand = $self->ligand();

        my $key_atoms1 = shift @key_atoms1;
        my $key_atoms2 = shift @key_atoms2;

        if (@$ref1_center == 2) {
            my ($subs, $ligand_mapped) =
                $ligand->_map_remote($ligand_ref, $key_atoms1, $key_atoms2, $ref1_center);
            $self->replace_ligand($ligand_mapped);
            $self->minimize_sub_torsions(object => $self->ligand(), subs => $subs);
        }

        while (@key_atoms1) {
            my $ligand = $self->ligand();
            my $key_atoms1 = shift @key_atoms1;
            my $key_atoms2 = shift @key_atoms2;

            my ($subs, $ligand_mapped) =
                $ligand->_map_remote($ligand_ref, $key_atoms1, $key_atoms2, $ref1_center);
            $self->minimize_sub_torsions(object => $self->ligand(), subs => $subs);
        }

        $self->ligand()->RMSD( ref_geo => $ligand_ref,
                               ref_atoms1 => [@key_total1],
                               ref_atoms2 => [@key_total2] );
    }else {
        $self->minimize_sub_torsions(object => $self->ligand());
    }
    $self->rebuild();
}


sub _map_center {
    my $self = shift;

    my ($ligand, $ref1, $ref2) = @_;
    my $ligand_ref = $self->ligand();

    #first block, this is a little more complex
    my $single_center;
    if (@$ref1 == 1) {
        $single_center = 1;

        ($ref1, $ref2) = $ligand->_explore_ref_atoms($ligand_ref,
                                                     $ref1,
                                                     $ref2,
                                                     1);
    }

    if (@$ref1 > 3 && $single_center) {
        my ($noC1) = grep{ $ref1->[$_] ne 'C' } (0..$#{$ref1});
        my ($noC2) = grep{ $ref2->[$_] ne 'C' } (0..$#{$ref2});

        $noC1 //= 0;
        $noC2 //= 0;

        my @ref1_temp = @$ref1[$noC1..$#{$ref1}];
        my @ref2_temp = @$ref2[$noC2..$#{$ref2}];

        push (@ref1_temp, @$ref1[1..$noC1-1]);
        push (@ref2_temp, @$ref2[1..$noC2-1]);

        unshift(@ref1_temp, $ref1->[0]);
        unshift(@ref2_temp, $ref2->[0]);


        unless ($noC1 && $noC2) {
            #try mapping multiple times
            my $E_min = 1E20;
            my $ligand_min;
            my $active_center = shift @$ref1;
            for my $i (0..$#{$ref1}) {
                my @ref1_temp = @$ref1[$i..$#{$ref1}];
                push(@ref1_temp, @$ref1[0..$i-1]);
                unshift(@ref1_temp, $active_center);
                my $ligand_copy = $ligand->copy();
                $ligand_copy->RMSD(ref_geo => $ligand_ref,
                                   ref_atoms1 => [@ref1_temp],
                                   ref_atoms2 => [@ref2_temp]);

                $self->replace_ligand($ligand_copy);

                my $LJ = $self->LJ_energy($self->ligand());
                if ($LJ < $E_min) {
                    $E_min = $LJ;
                    $ligand_min = $ligand_copy;
                }
            }

            $self->replace_ligand($ligand_min);
        }else {
            $ligand->RMSD( ref_geo => $ligand_ref,
                           ref_atoms1 => [@ref1_temp],
                           ref_atoms2 => [@ref2_temp] );
            $self->replace_ligand($ligand);
        }

        return [@ref1_temp];
    }elsif (@$ref1 > 2) {
        $ligand->RMSD( ref_geo => $ligand_ref,
                       ref_atoms1 => $ref1,
                       ref_atoms2 => $ref2 );
        $self->replace_ligand($ligand);

        return $ref1;
    }elsif (@$ref1 == 2 &&
            (@$ref2 == 2)) {
        $ligand->RMSD( ref_geo => $ligand_ref,
                       ref_atoms1 => $ref1,
                       ref_atoms2 => $ref2 );
        $self->replace_ligand($ligand);

        my $point = $self->ligand()->get_point($ref1->[0]);
        my $axis;
        if ($ref1->[1] =~ /^\d+$/) {
            $axis = $self->ligand()->get_bond($ref1->[1], $ref1->[0]);
        }else {
            $axis = $ref1->[1] - $self->ligand()->get_point($ref1->[0]);
        }

        $self->minimize_rotation( object => $LIGAND,
                                  point => $point,
                                  axis => $axis );
        return $ref1;
    }
}


#This method screen a virtual libary of potential catalyst design
sub screen_subs {
    my $self = shift;

    my $component = shift;

    my %params = @_;

    my $object;
    if ($component eq $LIGAND) { $object = $self->ligand() };
    if ($component eq $SUBSTRATE) { $object = $self->substrate() };

	# ligand can take named substitutions, substrate cannot
	my @subs;
	if ( $component eq $LIGAND ) {
		@subs = keys %params;
	} else {
		@subs = grep { $_ =~ /(\d+,?)+/ } keys %params;
	}

    my @subs_final = ({});
    for my $key (@subs) {
        my @subs_temp = map { {%{ $_ }} } @subs_final;
        @subs_final = ();
        for my $sub (@{ $params{$key} }) {
            my @subs_temp_temp = map { {%{ $_ }} } @subs_temp;
            map { my $sub_temp = $_;
                  map { $sub_temp->{$_} = $sub } split(/,/, $key) }
                @subs_temp_temp;
            push (@subs_final, @subs_temp_temp);
        }
    }

    my @cata = map{$self->copy()} (0..$#subs_final);

    map {$cata[$_]->substitute($component, %{ $subs_final[$_] })} (0..$#subs_final);

    return @cata;
}


#This is a substitute method for the catalysis class,
#you just provide atom numbers of the catalysis system
#and the substitutent you want to make.
#This can only be used when your substitutent target is also a
sub substitute {
    my $self = shift;

    my $component = shift;

    my %substituents = @_;

    $component //= '';

    my $object;
    if ($component eq $LIGAND) { $object = $self->ligand() };
    if ($component eq $SUBSTRATE) { $object = $self->substrate() };

    if (! $object) {
        $self->_replace_all();
        $self->remove_clash;
        $self->rebuild();
        return;
    }

    my %subs_final;

    my @inexplicit_subs = grep { $_ !~ /^\d+$/ } keys %substituents;
    my @explicit_subs = grep { $_ =~ /^\d+$/ } keys %substituents;

    for my $inexplicit_sub (@inexplicit_subs) {
        my @subs = grep { $object->{substituents}->{$_}->{name} eq $inexplicit_sub }
                        keys %{ $object->{substituents} };
        map { $subs_final{$_} = $substituents{$inexplicit_sub} } @subs;
    }

    map { $subs_final{$_} = $substituents{$_} } @explicit_subs;

    $object->_substitute(%subs_final);

    $self->_replace_all();

    $self->remove_clash;
    $self->rebuild();
}


sub _replace_all {
    my $self = shift;
    for my $object ($self->ligand(), $self->substrate()) {
        for my $target ( sort { $a <=> $b } keys %{ $object->{substituents} } ) {
            next unless $object->{substituents}->{$target}->{sub};
            next if ($object->{substituents}->{$target}->{sub} =~ /^\d+(,\d+)*$/);
            $object->replace_substituent(
                target => $target,
                   sub => $object->{substituents}->{$target}->{sub} );
            $self->rebuild_coords();
            $self->minimize_sub_torsion( object => $object,
                                         target => $target );
            delete $object->{substituents}->{$target}->{sub};
            #mark as the user defined;
        }
    }
}


sub replace_ligand {
    my $self = shift;

    my ($ligand) = @_;

    my $self_l_copy = $self->ligand()->copy();

    #replace the atom in the constraints with the key atoms in the new ligand
    #first put the active centers in a one dimensional array;
    my @self_l_key = map { @$_ } @{ $self_l_copy->{active_centers} };
    my @ligand_key = map { @$_ } @{ $ligand->{active_centers} };
    for my $con_n (0..$#{ $self_l_copy->{constraints} }) {
        for my $i (0..1) {
            for my $center_n (0..$#self_l_key) {
                if ($self_l_copy->{constraints}->[$con_n]->[$i] =~ /^d+$/ &&
                    $self_l_copy->{constraints}->[$con_n]->[$i] == $self_l_key[$i]) {
                    $self_l_copy->{constraints}->[$con_n]->[$i] = $ligand_key[$i];
                }
            }
        }
    }

    $ligand->{constraints} = $self_l_copy->{constraints};

    $self->{ligand} = $ligand;
    $self->rebuild();
}



#rotate one part of the catalysis to minimize the RMSD of the whole system
#point is a ->point return
#v is a vector of the axis
sub minimize_rotation {
    my $self = shift;
    my %param = @_;
    my ($point, $axis) = ($param{point}, $param{axis});

    my $object;
    if ($param{object} eq $LIGAND) {$object = $self->ligand()};
    if ($param{object} eq $SUBSTRATE) {$object = $self->substrate()};

    my $increment = 5;

    my $E_min = 1E20;
    my $angle_min = 0;

    foreach my $count (0..360/$increment + 1) {
        my $angle = $count*$increment;
        $object->center_genrotate($point, $axis, deg2rad($increment));
        my $energy = $self->LJ_energy($object);
        if($energy < $E_min) {
          $angle_min = $angle;
          $E_min = $energy;
        }
    }

    $object->center_genrotate($point, $axis, deg2rad($angle_min));
    $self->rebuild_coords();
}


sub minimize_sub_torsions {
    my $self = shift;
    my %params = @_;
    my ($object, $subs) = ($params{object}, $params{subs});

    $subs //= [ keys %{ $object->{substituents} } ];

    for my $target ( @$subs ) {
        $self->minimize_sub_torsion(object => $object, target => $target);
    }
}


#rotate one part of the catalysis to minimize the RMSD of the whole system
#point is a ->point return
#v is a vector of the axis
sub minimize_sub_torsion {
    my $self = shift;
    my %param = @_;
    my ($object, $target) = ($param{object}, $param{target});

    my $increment = 15;

    my $E_min = 1E20;
    my $angle_min = 0;

    foreach my $count (0..360/$increment + 1) {
        my $angle = $count*$increment;
        $object->sub_rotate( target => $target,
                              angle => $increment );
        my $energy = $self->part_LJ_energy($object->{substituents}->{$target});
        if($energy < $E_min) {
          $angle_min = $angle;
          $E_min = $energy;
        }
    }

    $object->sub_rotate( target => $target, angle => $angle_min);
}


sub part_LJ_energy {
    my $self = shift;

    my $part = shift;

    my @rest;

    for my $i (0..$#{$self->{elements}}) {
        unless (grep {@{$self->{coords}->[$i]} ~~ @$_} @{$part->{coords}}) {
            push (@rest, $i);
        }
    }

    my $rest_geo = $self->subgeo([@rest]);

    return $rest_geo->LJ_energy_with($part);
}


#return on which component the atom is
sub on_which_component {
    my ($self, $atom) = @_;

    my $return;
    if (grep { $_ == $atom } @{ $self->{substrate_atoms} }) {
        $return = $self->substrate();
    }
    if ($return == $self->{center_atom}) {
        $return = $self->center();
    }
    if (grep { $_ == $atom } @{ $self->{ligand_atoms} }) {
        $return = $self->ligand();
    }

    return $return;
}


sub number_of_conformers {
    my $self = shift;
    my %params = @_;
    my ($object) = ( $params{component} );

    $object //= '';

    my $cf_number_substrate = 1;
    my $cf_number_ligand = 1;

    for my $sub (keys %{ $self->substrate()->{substituents} }) {
        $cf_number_substrate *= $self->substrate()->{substituents}->{$sub}->{conformer_num};
    }

    for my $sub (keys %{ $self->ligand()->{substituents} }) {
        $cf_number_ligand *= $self->ligand()->{substituents}->{$sub}->{conformer_num};
    }

    if ($object eq $SUBSTRATE) {
        return $cf_number_substrate;
    }elsif ($object eq $LIGAND) {
        return $cf_number_ligand;
    }else {
        my $total = $cf_number_substrate * $cf_number_ligand;
        return $total;
    }
}


#FIXME make conformer for the whole geometry
#not only ligand or substratwe
sub make_conformer {
    my $self = shift;
    my %params = @_;
    my ($number, $current_number,
        $new_array, $old_array) = ( $params{new_number},
                                    $params{current_number},
                                    $params{new_array},
                                    $params{current_array} );

    $current_number //= $self->{conf_num};

    $new_array //= $self->conf_array(number => $number);
    $old_array //= $self->conf_array(number => $current_number);

    my $num_l = keys %{ $self->ligand()->{substituents} };
    my $num_s = keys %{ $self->substrate()->{substituents} };

    if ($num_l > 0) {
        $self->ligand()->make_conformer( current => [@$old_array[0..$num_l - 1]],
                                                new => [@$new_array[0..$num_l - 1]] );
    }
    if ($num_s > 0) {
        $self->substrate()->make_conformer( current => [@$old_array[$num_l..$num_l + $num_s - 1]],
                                             new => [@$new_array[$num_l..$num_l + $num_s - 1]] );
    }
    $self->rebuild();
    $self->{conf_num} = $number;
}


sub conf_array {
    my $self = shift;

    my %params = @_;

    my ($number) = ($params{number});

    $number //= $self->{conf_num};

    my $conf_max_array = $self->_conf_max_array();

    my @reverse_max = reverse @$conf_max_array;

    my @conf_array;

    for my $max (@reverse_max) {
        my $temp = (($number - 1) % $max) + 1;
        push (@conf_array, $temp);
        $number = int(($number - 1)/$max) + 1;
    }

    return ([reverse @conf_array]);
}


sub conf_num {

    my $self = shift;

    my ($array) = @_;

    my @array = @$array;

    my $conf_max_array = $self->_conf_max_array();

    my $conf_num = $array[0];

    for my $i (1..$#array) {
        $conf_num = ($conf_num - 1) * $conf_max_array->[$i] + $array[$i];
    }

    return $conf_num;
}


sub _conf_max_array {
    my $self = shift;

    my @conf_max_array;

    for my $target (sort { $a <=> $b } keys %{ $self->ligand()->{substituents} }) {
        #NOTE reverse order here
        push(@conf_max_array, $self->ligand()->{substituents}->{$target}->{conformer_num});
    }

    for my $target (sort { $a <=> $b } keys %{ $self->substrate()->{substituents} }) {
        #NOTE reverse order here
        push(@conf_max_array, $self->substrate()->{substituents}->{$target}->{conformer_num});
    }

    return [@conf_max_array];
}


sub max_conf_number {
    my $self = shift;

    my ($sub) = @_;
    my $max;

    my @substrate_subs = sort { $a <=> $b } keys %{ $self->substrate()->{substituents} };
    my @ligand_subs = sort { $a <=> $b } keys %{ $self->ligand()->{substituents} };

    if ($sub <= $#ligand_subs) {
        $max = $self->ligand()->{substituents}->{$ligand_subs[$sub]}->{conformer_num};
    }else {
        $max = $self->substrate()->{substituents}->{$substrate_subs[$sub-$#ligand_subs-1]}->{conformer_num};
    }

    return $max;
}


sub init_conf_num {

    my $self = shift;

    $self->{conf_num} = 1;
}


sub conformer_rmsd {
    my $self = shift;

    my %params = @_;

    my ($cata2, $heavy_only) = ( $params{ref_cata}, $params{heavy_atoms} );

    $heavy_only //= 0;

    #get substituents

    my @sub_subs_keys = keys %{ $self->substrate()->{substituents} };
    my @lig_subs_keys = keys %{ $self->ligand()->{substituents} };

    my @sub_subs_atoms = map { $self->substrate()->{substituents}->{$_}->{coords} }
                            @sub_subs_keys;
    my @lig_subs_atoms = map { $self->ligand()->{substituents}->{$_}->{coords} }
                            @lig_subs_keys;

    my @subs_atoms_nums;
    for my $sub_atoms (@sub_subs_atoms, @lig_subs_atoms) {
        my @subs_nums;
        for my $atom(0..$#{$sub_atoms}) {
            my ($atom_num) = grep { $self->{coords}->[$_] == $sub_atoms->[$atom] }
                                (0..$#{ $self->{coords}});
            push (@subs_nums, $atom_num);
        }
        push (@subs_atoms_nums, [@subs_nums]);
    }

    my $cata1 = $self->copy();

    my @rmsd;
    for my $reorder (0,1) {
        if ($reorder) {
            for my $sub_atoms_nums (@subs_atoms_nums) {
                my %visited;
                while (@$sub_atoms_nums) {
                    for (my $i=0; $i<=$#{$sub_atoms_nums}; $i++) {
                        next if (exists $visited{$sub_atoms_nums->[$i]});

                        my $short_d = 999;
                        my $atom_i = $sub_atoms_nums->[$i];
                        my $short_atom = $atom_i;
                        my $short_num = $i;
                        for my $j (0..$#{$sub_atoms_nums}) {
                            #next if (exists $visited{$sub_atoms_nums->[$j]});
                            my $atom_j = $sub_atoms_nums->[$j];
                            if($cata1->{elements}->[$atom_i] eq $cata2->{elements}->[$atom_j]) {
                                my $distance = $cata2->distance( atom1 => $atom_i,
                                                                 atom2 => $atom_j,
                                                             geometry2 => $cata1 );
                                if ($distance < $short_d) {
                                    $short_d = $distance;
                                    $short_atom = $atom_j;
                                    $short_num = $j;
                                }
                            }
                        }
                        $visited{$atom_i} = 1;

                        next if ($atom_i == $short_atom);
                        if ($short_atom < $atom_i) {
                            my $current_d = $cata2->distance( atom1 => $atom_i,
                                                              atom2 => $atom_i,
                                                          geometry2 => $cata1 )
                                            +
                                            $cata2->distance( atom1 => $short_atom,
                                                              atom2 => $short_atom,
                                                          geometry2 => $cata1 );

                            my $swap_d = $cata2->distance( atom1 => $atom_i,
                                                           atom2 => $short_atom,
                                                       geometry2 => $cata1 )
                                         +
                                         $cata2->distance( atom1 => $short_atom,
                                                           atom2 => $atom_i,
                                                       geometry2 => $cata1 );

                            next if ($current_d < $swap_d);
                        }

                        my $temp = $cata1->{coords}->[$atom_i];
                        $cata1->{coords}->[$atom_i] = $cata1->{coords}->[$short_atom];
                        $cata1->{coords}->[$short_atom] = $temp;
                        $i = $short_num - 1;
                    }

                    $sub_atoms_nums = [ grep {!exists $visited{$_}} @$sub_atoms_nums ];
                }
            }
        }

        push (@rmsd, $cata1->RMSD(ref_geo => $cata2));
    }

    @rmsd = sort { $a <=> $b } @rmsd;

    return $rmsd[0];
}


sub remove_clash {
    my $self = shift;
    my $CRASH_CUTOFF = $CUTOFF->{CRASH_CUTOFF};

    my $relief = 1;

    for my $object ($self->substrate(), $self->ligand()) {
        for my $key (keys %{ $object->{substituents} }) {
            my %crash;
            my $axis;
            my $sub = $object->{substituents}->{$key};

            my $get_crash = sub {
                %crash = ();
                for my $atom_sub (0.. $#{ $object->{substituents}->{$key}->{elements} }) {
                    my @crash = grep {$self->distance( atom1 => $_,
                                                       atom2 => $atom_sub,
                                                   geometry2 => $sub)
                                            < $CRASH_CUTOFF} (0..$#{ $self->{elements} });

                    my ($atom_num) = grep { $self->{coords}->[$_] == $sub->{coords}->[$atom_sub] }
                                        (0..$#{ $self->{coords}});
                    @crash = grep {$_ != $atom_num} @crash;

                    @crash{@crash} = ();
                }

                my $vector_sum = V(0,0,0);
                map {$vector_sum += $self->get_bond($_, $sub->{end}, $object)} keys %crash;
                my $bond = $object->get_bond($key, $sub->{end});
                $axis = $vector_sum x $bond;
            };

            &$get_crash();

            my $b_times = 0;
            my @r_angles = (5, -10, 15, -20, 25, -30);
            my $end_point = $object->get_point($sub->{end});
            while(%crash) {
                if ($b_times < 3) {
                    $sub->center_genrotate($end_point, $axis, deg2rad(5));
                    $b_times ++;
                }elsif (@r_angles) {
                    $object->sub_rotate( target => $key, angle => shift @r_angles);
                }else {
                    $relief = 0;
                    last;
                }
                &$get_crash();
            }
        }
    }
    return $relief;
}


package AaronTools::Component;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
use Data::Dumper;
our @ISA = qw(AaronTools::Geometry);

sub new {
    my $class = shift;
    my %params = @_;

    my $substituents = $params{substituents};

    my $self = new AaronTools::Geometry();
    bless $self, $class;

    $self->{substituents} = {};

    if (exists $params{name}) {
        $self->set_name($params{name});
    }

    $self->set_substituents($substituents);

    return $self;
}


sub backbone {
    my $self = shift;
    return $self->{backbone};
}


sub get_all_sub_targets {
    my $self = shift;
    my @targets = sort { $a <=> $b } keys %{ $self->{substituents} };
    return [@targets];
}


#return the atom number of the atom where the substituent is on
#the backbone
sub get_substituent {
    my ($self, $target) = @_;
    return $self->{substituents}->{$target};
}


sub get_substituent_atoms {
    my $self = shift;
    my ($target) = @_;
    my $atom2 = $self->get_substituent($target)->end();
    return $self->get_all_connected($target, $atom2);
}


#set substituents from a substiuents information hash
sub set_substituents {

    my ($self, $substituents) = @_;

    if ($substituents) {
        for my $key (keys %{ $substituents }) {
            $self->{substituents}->{$key} = new AaronTools::Substituent(
                                                sub => $substituents->{$key} );
            if ($self->{backbone}) {
                my @connected = @{$self->{connection}->[$key]};
                my @back = grep { $_ <= $#{ $self->{backbone}->{elements} } } @connected;
                my ($end) = @back
                                ? @back == 1
                                    ? @back
                                    : $self->get_sub($key)
                                : ('');

                $self->{substituents}->{$key}->{end} = $end;
            }
        }
    }
}


sub sub_rotate {
    my $self = shift;

    my %params = @_;

    my ($target, $angle) = ( $params{target}, $params{angle} );

    my $sub = $self->{substituents}->{$target};

    my $point = $self->get_point($sub->{end});
    my $axis = $self->get_bond($sub->{end}, $target);

    $sub->center_genrotate($point, $axis, deg2rad($angle));
}


#rebuid self geometry from backbone and substituents
sub rebuild {
    my $self = shift;

    my @elements = @{ $self->{backbone}->{elements} };
    my @flags = @{ $self->{backbone}->{flags} };
    my @coords = @{ $self->{backbone}->{coords} };

    if ($self->{substituents}) {
        for my $target (sort {$a <=> $b} keys %{ $self->{substituents} }) {
            my @element = @{ $self->{substituents}->{$target}->{elements} };
            $self->{backbone}->{elements}->[$target] = $elements[$target] = shift @element;
            my @flag = @{ $self->{substituents}->{$target}->{flags} };
            $self->{backbone}->{flags}->[$target] = $flags[$target] = shift @flag;
            my @coord = @{ $self->{substituents}->{$target}->{coords} };
            $self->{backbone}->{coords}->[$target] = $coords[$target] = shift @coord;
            push (@elements, @element);
            push (@flags, @flag);
            push (@coords, @coord);
        }
    }

    $self->{elements} = [@elements];
    $self->{flags} = [@flags];
    $self->{coords} = [@coords];

    $self->refresh_connected();
}


sub freeze_all_atoms {
    my $self = shift;

    $self->_set_flags_all_atoms(-1);
}


sub relax_all_atoms {
    my $self = shift;

    $self->_set_flags_all_atoms(0);
}


sub _set_flags_all_atoms {
    my $self = shift;
    my ($flag) = @_;

    $self->backbone()->{flags} = [ map { $flag } @{ $self->backbone()->{flags} } ];

    for my $target (keys %{ $self->{substituents} }) {
        $self->{substituents}->{$target}->{flags} = [ map { $flag }
                                    @{ $self->{substituents}->{$target}->{flags} } ]
    }

    $self->{flags} = [ map { $flag } @{ $self->{flags} } ];
}


#This is to update the backbone and substituents from the ligand or substrate coords,
#NOTE This is only desinged for different geometry with only internal manipulation
#The end of this function is to garantee the coords and substituents are the same array
#ref as in the coords
sub update_geometry {
    my $self = shift;

    $self->{backbone}->{coords} =
        [ @{ $self->{coords} }[0..$#{ $self->{backbone}->{coords} }] ];

    my $head = $#{ $self->{backbone}->{coords} } + 1;
    for my $target (sort {$a <=> $b} keys %{ $self->{substituents} }) {
        my $coords = [$self->{coords}->[$target],
                      @{ $self->{coords} }[$head..$head +
                        $#{ $self->{substituents}->{$target}->{elements} } - 1]];
        $self->{substituents}->{$target}->{coords} = $coords;
        $head += $#{ $coords };
    }
}


sub copy {
    my $self = shift;

    my $new =  new AaronTools::Geometry();

    bless $new, "AaronTools::Component";

    $new->{name} = $self->{name};
    $new->{constraint} = [ map { [ @$_ ] } @{ $self->{constraint} } ];

    if ($self->{substituents}) {
        $new->{substituents} = {};
        for my $target ( sort { $a <=> $b } keys %{ $self->{substituents} } ) {
            $new->{substituents}->{$target} =
                $self->{substituents}->{$target}->copy();
        }
    }

    if ($self->{backbone}) {
        $new->{backbone} = $self->{backbone}->copy();
    }

    if ($self->{constraints}) {
        $new->{constraints} = [ map { [ @$_ ] } @{ $self->{constraints} }];
    }

    if ($self->{backbone}) {
        $new->rebuild();
    }else {
        $new->{elements} = [ @{ $self->{elements} } ];
        $new->{flags} = [ @{ $self->{flags} } ];
        $new->{coords} = [ map { [ @$_ ] } @{ $self->{coords} }];
    }
    $new->refresh_connected();

    return $new;
}


#This is to replace an old substituent
sub replace_substituent {
    my $self = shift;

    my %params = @_;
    my ($sub, $target) = ($params{sub}, $params{target});
    my $end = $self->{substituents}->{$target}->{end};

    my $new_sub = new AaronTools::Substituent( name => $sub, end => $end );

    $new_sub->_align_on_geometry( geo => $self, target => $target, end => $end );

    $new_sub->{end} = $self->{substituents}->{$target}->{end};

    $new_sub->{flags} = [ map { 0 } @{ $new_sub->{flags} } ];

    $self->{substituents}->{$target} = $new_sub;
}

sub _rearrange_active_con {
    my ($self, $old_new) = @_;

    if ($self->{constraints}) {
        for my $n_con (0..$#{ $self->{constraints} }) {
            for my $i (0..1) {
                if ($self->{constraints}->[$n_con]->[$i] =~ /^\d+$/ &&
                    exists $old_new->{$self->{constraints}->[$n_con]->[$i]}) {
                    $self->{constraints}->[$n_con]->[$i] =
                        $old_new->{$self->{constraints}->[$n_con]->[$i]};
                }
            }
        }
    }

    if ($self->{substituents}) {
        for my $key (sort {$a <=> $b} keys %{ $self->{substituents} }) {
            if (exists $old_new->{$self->{substituents}->{$key}->{end}} &&
                ($old_new->{$self->{substituents}->{$key}->{end}} !=
                $self->{substituents}->{$key}->{end})) {
                $self->{substituents}->{$key}->{end} = $old_new->{$self->{substituents}->{$key}->{end}};
            }
            if (exists $old_new->{$key} &&
                ($old_new->{$key} != $key)) {
                $self->{substituents}->{$old_new->{$key}} = $self->{substituents}->{$key};
                delete $self->{substituents}->{$key};
            }
        }
    }

    if ($self->{active_centers}) {
        for my $i (0..$#{ $self->{active_centers} }) {
            $self->{active_centers}->[$i] =
                [ map { exists $old_new->{$_} ? $old_new->{$_} : $_ }
                      @{ $self->{active_centers}->[$i] } ];
        }
    }
}


sub make_conformer {
    my $self = shift;

    my %params = @_;

    my ($new_array, $old_array) = ( $params{new}, $params{current});
    my @targets = sort { $a <=> $b } keys %{ $self->{substituents} };

    for my $i (0..$#targets) {
        my $angle = $self->{substituents}->{$targets[$i]}->{conformer_angle} *
                    ( $new_array->[$i] - $old_array->[$i] );
        $self->sub_rotate( target => $targets[$i], angle => $angle );
        $self->{substituents}->{$targets[$i]}->{flags} = [(0) x
            scalar @{ $self->{substituents}->{$targets[$i]}->{flags} }];
    }
}


sub _detect_substituent {
    my $self = shift;
    my %params = @_;

    my ($start, $end) = ($params{target}, $params{end});

    my $is_sub = $self->__detect_substituent($start, $end);

    if ($is_sub) {
        my $to_sub = $self->{substituents}->{$start}->{sub};
        $self->{substituents}->{$start} = $self->detect_substituent( target => $start,
                                                                        end => $end );
        $self->{substituents}->{$start}->{sub} = $to_sub;
    }

    return $is_sub;
}


sub __detect_substituent {
    my $self = shift;

    my ($start, $end) = @_;

    my $connected = $self->get_all_connected($start, $end);

    my $is_sub = 1;

    for my $sub_target (sort {$a <=> $b} keys %{ $self->{substituents} }) {
        if ($start != $sub_target &&
            grep { $_ == $sub_target } @$connected) {
            $is_sub = 0;

            #set substituents
            my $to_sub = $self->{substituents}->{$sub_target}->{sub};
            my ($end, $old_sub_atoms) = $self->get_sub($sub_target);
            $self->{substituents}->{$sub_target} = $self->detect_substituent(
                                                        target => $sub_target,
                                                           end => $end );
            $self->{substituents}->{$sub_target}->{sub} = $to_sub;

            last;
        }
    }

    return $is_sub;
}


#This is the not-called get all connected, it will return
#atoms connected to start atoms but exclude all substituents in the
#substituents attribute
sub _get_all_connected {
    my ($self, $start, $end) = @_;

    my $connected = $self->get_all_connected($start, $end);

    my %connected_hash;
    @connected_hash{ @$connected } = ();

    my @sub_keys = sort {$a <=> $b} keys %{ $self->{substituents} };

    @sub_keys = grep { $_ != $start } @sub_keys;

    for my $sub_target (@sub_keys) {
        if (grep { $_ == $sub_target } @$connected) {
            my $connected_sub = $self->get_all_connected(
                                    $sub_target,
                                    $self->{substituents}->{$sub_target}->{end} );
            my @connected_sub = @$connected_sub[1..$#{ $connected_sub }];

            delete @connected_hash{ @connected_sub };
        }
    }
    my $connected_atoms = [ sort { $a <=> $b } keys %connected_hash ];

    return ($connected_atoms);
}


sub substitute {
    my $self = shift;

    $self->_substitute(@_);

    $self->_replace_all();
}


sub _substitute {
    my $self = shift;
    my %substituents = @_;

    my %new_subs;
    my %old_subs;
    my @current_subs = keys %{ $self->{substituents} };
    my @new_subs = keys %substituents;
    @new_subs{@new_subs} = ();
    @old_subs{@new_subs} = ();
    delete @new_subs{@current_subs};
    delete @old_subs{keys %new_subs};

    map { $self->{substituents}->{$_}->{sub} = $substituents{$_} } keys %old_subs;

    unless (%new_subs) {
        return;
    }

    %substituents = map { $_ => $substituents{$_} } keys %new_subs;

    $self->set_substituents(\%substituents);

    my $old_new = {};
    my $refresh_backbone;

    my @subs_with_end = grep {$self->{substituents}->{$_}->{end}}
                            keys %{$self->{substituents}};

    for my $sub_key (sort { $a <=> $b } @subs_with_end) {
        my $end = $self->{substituents}->{$sub_key}->{end};
        if (! $self->__detect_substituent($sub_key, $end)) {
            my $connected_atoms = $self->_get_all_connected($sub_key, $end);

            my $back_atoms = $#{ $self->{backbone}->{elements} } + 1;

            map { $old_new->{$connected_atoms->[$_]} = $back_atoms + $_ - 1 }
                (1..$#{$connected_atoms});

            for my $i (1..$#{ $connected_atoms }) {
                for my $att ('elements', 'flags', 'coords') {
                    my $temp = $i+$back_atoms;
                    $self->{backbone}->{$att}->[$i + $back_atoms - 1] =
                        $self->{$att}->[$connected_atoms->[$i]];
                }
            }

            delete $self->{substituents}->{$sub_key};
        }elsif ( ! @{ $self->{substituents}->{$sub_key}->{elements} } ) {
            my $to_sub = $self->{substituents}->{$sub_key}->{sub};
            $self->{substituents}->{$sub_key} = $self->detect_substituent(
                                                        target => $sub_key,
                                                           end => $end );
            $self->{substituents}->{$sub_key}->{sub} = $to_sub;

            for my $coord (@{ $self->{substituents}->{$sub_key}->{coords} }) {
                $refresh_backbone = grep { @$_ ~~ @$coord }
                                        @{ $self->{backbone}->{coords} };
                last if $refresh_backbone;
            }
        }
    }
    $self->_rearrange_active_con($old_new);

    $self->detect_backbone_subs() if $refresh_backbone;

    $self->rebuild();
}


sub _replace_all {
    my $self = shift;

    my ($minimize) = @_;

    for my $target ( sort { $a <=> $b } keys %{ $self->{substituents} } ) {
        if ($self->{substituents}->{$target}->{sub}) {
            $self->replace_substituent(
                target => $target,
                   sub => $self->{substituents}->{$target}->{sub} );
            $self->rebuild();
            if ($minimize) {
                $self->minimize_sub_torsion( target => $target );
            }
            delete $self->{substituents}->{$target}->{sub};
            #mark as the user defined;
        }
    }
}


#rotate one part of the catalysis to minimize the RMSD of the whole system
#point is a ->point return
#v is a vector of the axis
sub minimize_sub_torsion {
    my $self = shift;
    my ($target) = @_;

    my $increment = 5;

    my $E_min = 1E20;
    my $angle_min = 0;

    foreach my $count (0..360/$increment + 1) {
        my $angle = $count*$increment;
        $self->sub_rotate( target => $target,
                            angle => $increment );
        my $energy = $self->part_LJ_energy($self->{substituents}->{$target});
        if($energy < $E_min) {
          $angle_min = $angle;
          $E_min = $energy;
        }
    }

    $self->sub_rotate( target => $target, angle => $angle_min);
}


sub part_LJ_energy {
    my $self = shift;

    my $part = shift;

    my @rest;

    for my $i (0..$#{$self->{elements}}) {
        unless (grep {@{$self->{coords}->[$i]} ~~ @$_} @{$part->{coords}}) {
            push (@rest, $i);
        }
    }

    my $rest_geo = $self->subgeo([@rest]);

    return $rest_geo->LJ_energy_with($part);
}



package AaronTools::Ligand;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
use Data::Dumper;
our @ISA = qw(AaronTools::Component AaronTools::Geometry);

sub new {
    my $class = shift;
    my %params = @_;

    my $self = new AaronTools::Component( substituents => $params{substituents} );
    bless $self, $class;

    if (exists $params{name}) {
        $self->set_name($params{name});
        if (-f "$QCHASM/AaronTools/Ligands/$self->{name}.xyz") {
            $self->read_geometry("$QCHASM/AaronTools/Ligands/$self->{name}.xyz");
        }

    }

    if ($self->{active_centers}) {
        $self->detect_backbone_subs( no_new_subs => $params{no_new_subs} );
    }

    return $self;
}


#This function detect backbone and substituents connected on the backbone
#$self->{substituents} and $self->{backbone} will be generated.
sub detect_backbone_subs {
    my $self = shift;

    my %params = @_;

    my $no_new_subs = $params{no_new_subs};

    my @active_centers = map {@$_} @{ $self->{active_centers} };

    my %backbone = %{ $self->bare_backbone([@active_centers]) };

    for my $atom (keys %backbone) {
        for my $atom_connected (@{ $self->{connection}->[$atom] }) {
            if (!exists $backbone{$atom_connected}) {
                $backbone{$atom_connected} = 1;
                #if double or triple connected to the backbone
                unless ($#{ $self->{connection}->[$atom_connected] } == 0 &&
                        ($#{ $self->{connection}->[$atom_connected] } + 1 <
                        $CONNECTIVITY->{$self->{elements}->[$atom_connected]})) {
                    if ($self->{elements}->[$atom_connected] ne 'H') {
                        #dectec if the substituent contains any substituents has been
                        #made before. If yes, including this substituent in the
                        #backbone, if not make a new backbone
                        my $is_sub = 1;
                        if ($no_new_subs) {
                            unless (grep { $_ == $atom_connected}
                                        keys %{ $self->{substituents} }) {
                                $is_sub = 0;
                            }
                        }

                        unless ($is_sub &&
                                $self->_detect_substituent(target => $atom_connected,
                                                           end => $atom)) {
                            my $new_back_atoms = $self->_get_all_connected( $atom_connected,
                                                                            $atom );
                            for my $atom_temp (@$new_back_atoms) {
                                $backbone{$atom_temp} = 1;
                            }
                        }
                    }elsif ($self->{substituents}->{$atom_connected}) {
                            $self->_detect_substituent(target => $atom_connected,
                                                       end => $atom);
                    }
                }
            }
        }
    }

    my @final_backbone = sort { $a <=> $b } keys %backbone;

    $self->{backbone} = $self->subgeo([@final_backbone]);
    #rearrange active centers and constraints
    my $backbone_capped = { map { $final_backbone[$_] => $_ } (0..$#final_backbone) };
    $self->_rearrange_active_con($backbone_capped);
    $self->rebuild();
}


sub capped_backbone {
    my $self = shift;

    my $geo = $self->backbone()->copy();
    for my $atom (keys %{ $self->{substituents} }) {
        $geo->{elements}->[$atom] = 'H';

        $geo->correct_bond_length( atom1 => $self->{substituents}->{$atom}->end(),
                                   atom2 => $atom );
    }
    return $geo;
}


sub copy {
    my $self = shift;

    my $new =  new AaronTools::Geometry( name => $self->{name},
                                         constraint => [ map { [ @$_ ] }
                                                            @{ $self->{constraint} } ] );
    bless $new, "AaronTools::Ligand";

    if ($self->{active_centers}) {
        $new->{active_centers} = [map { [@$_] } @{ $self->{active_centers} }];
    }

    $new->{RMSD_bonds} = [];
    if ($self->{RMSD_bonds}) {
        for my $i (0..$#{ $self->{RMSD_bonds} }) {
            $new->{RMSD_bonds}->[$i] = [ map { [@$_] } @{ $self->{RMSD_bonds}->[$i] }];
        }
    }

    if ($self->{backbone}) {
        $new->{backbone} = $self->{backbone}->copy();
    }

    if ($self->{substituents}) {
        $new->{substituents} = { map { $_ => $self->{substituents}->{$_}->copy() }
                                     keys %{ $self->{substituents} } };
    }

    if ($self->{constraints}) {
        $new->{constraints} = [ map { [ @$_ ] } @{ $self->{constraints} }];
    }
    $new->rebuild();
    $new->refresh_connected();

    return $new;
}


sub _explore_ref_atoms {
    my ($self, $ligand_ref, $ref_atoms1, $ref_atoms2, $first_block) = @_;

    my $atom1 = $ref_atoms1->[0];
    my $atom2 = $ref_atoms2->[0];

    my @connections1 = @{$self->{connection}->[$atom1]};
    my @connections2 = @{$ligand_ref->{connection}->[$atom2]};

    my $ref1 = [$atom1];
    my $ref2 = [$atom2];

    my @con1_back;
    my @con2_back;
    if (@connections2 > 1) {
        my $bare_back = $ligand_ref->bare_backbone($ref2);
        @con2_back = ();
        for my $atom (@connections2) {
            push (@con2_back, $atom) if (exists $bare_back->{$atom});
        }
    }

    if (@connections1 > 1) {
        my $bare_back = $self->bare_backbone($ref1);
        @con1_back = ();
        for my $atom (@connections1) {
            push (@con1_back, $atom) if (exists $bare_back->{$atom});
        }
    }

    @con1_back = @connections1 unless (@con1_back && (@connections1 == @connections2));
    @con2_back = @connections2 unless (@con2_back && (@connections1 == @connections2));

    if (@con1_back != @con2_back ||
        (@connections1 != @connections2) ||
        !$first_block) {
        my $center1 = $self->get_center([@con1_back]);
        my $center2 = $ligand_ref->get_center([@con2_back]);
        #Normalize the centers
        my $v1 = $center1 - $self->get_point($ref1->[0]);
        my $v2 = $center2 - $ligand_ref->get_point($ref2->[0]);

        $v1 *= 1/abs($v1);
        $v2 *= 1/abs($v2);

        $ref1->[1] = $center1 + $v1;
        $ref2->[1] = $center2 + $v2;

    }else {
        if (@con1_back == 2 &&
            (@con1_back ~~ @connections1)) {
            my $msg = "One active center in the backbone may cause incorrect mapping\n";
            warn($msg);

            push (@$ref1, @con1_back);
            push (@$ref2, @con2_back);
        }else {
            my $center1 = $self->get_center([@connections1]);
            my $center2 = $ligand_ref->get_center([@connections2]);
            my $axis1 = $center1 - $self->get_point($atom1);
            my $axis2 = $center2 - $ligand_ref->get_point($atom2);

            my $keep_sorting;

            do {
                $keep_sorting = 0;
                for my $i (0..$#con1_back-1) {
                    my $j = $i + 1;
                    my $vi = $self->get_point($con1_back[$i]) - $self->get_point($atom1);
                    my $vj = $self->get_point($con1_back[$j]) - $self->get_point($atom1);
                    if ($vi x $vj * $axis1 < 0) {
                        my $temp = $con1_back[$i];
                        $con1_back[$i] = $con1_back[$j];
                        $con1_back[$j] = $temp;
                        $keep_sorting = 1;
                    }
                }
            }while ($keep_sorting);

            do {
                $keep_sorting = 0;
                for my $i (0..$#con2_back-1) {
                    my $j = $i + 1;
                    my $vi = $ligand_ref->get_point($con2_back[$i]) -
                                $ligand_ref->get_point($atom2);
                    my $vj = $ligand_ref->get_point($con2_back[$j]) -
                                $ligand_ref->get_point($atom2);
                    if ($vi x $vj * $axis2 < 0) {
                        my $temp = $con1_back[$i];
                        $con2_back[$i] = $con1_back[$j];
                        $con2_back[$j] = $temp;
                        $keep_sorting = 1;
                    }
                }
            }while ($keep_sorting);

            push (@$ref1, @con1_back);
            push (@$ref2, @con2_back);
        }
    }

    return ($ref1, $ref2);
}


sub _map_remote {
    my $self = shift;

    my ($ref_ligand, $ref1, $ref2, $key_center) = @_;

    my $rotatable_bonds = $self->{backbone}->rotatable_bonds($ref1, $key_center);

    my @rotatable_bonds = map { [$_, $rotatable_bonds->{$_}] }
                            keys %{ $rotatable_bonds };

    my @fragments;
    for my $bond (@rotatable_bonds) {
        push (@fragments, $self->get_all_connected(@$bond));
    }

    my @subs;

    my $active = [@$key_center, @$ref1];

    my $backbone = $self->bare_backbone($active);

    if (@$ref1 == 1) {
        ($ref1, $ref2) = $self->_explore_ref_atoms($ref_ligand,
                                                   $ref1,
                                                   $ref2);
    }

    for my $sub (keys %{$self->{substituents}}) {
        if (grep {$_ == $self->{substituents}->{$sub}->{end}}
               keys %{$backbone}) {
            push (@subs, $sub);
        }
    }

    my $minRMSD = 9999;
    my $increment = 10;
    my $ligand;

    my $rmsd_align;

    my $last_rmsd = $minRMSD;
    my $keepgoing = 1;
    my $wrong_direction = 0;
    my $direction_changed = 0;
    my $angle = 0;

    my $rotate_angle = [(0) x @rotatable_bonds];
    my $min_rotate_angle;
    my $rotate_ith = -1;

    my $atoms_ref = [grep {$_ =~ /^\d+$/} @$ref1];

    $rmsd_align = sub {
        my ($bonds) = @_;
        $rotate_ith ++;
        my @bond = @{ shift @$bonds };
        my $point = $ligand->get_point($bond[1]);
        my $axis = $ligand->get_bond(reverse @bond);

        for my $n(1..360/$increment) {
            if (@$bonds) {
                my $bonds_temp = [ map { [@$_] } @$bonds ];
                &$rmsd_align($bonds_temp);
                $rotate_ith --;
            }
            $ligand->center_genrotate($point, $axis, deg2rad($increment), $atoms_ref);
            &__center_genrotate($ref1, $point, $axis, deg2rad($increment));
            $rotate_angle->[$rotate_ith] = $n*$increment;

            my $rmsd = $ligand->MSD( ref_geo => $ref_ligand,
                                   ref_atoms1 => $ref1,
                                   ref_atoms2 => $ref2,
                                   no_rot => 1);
            if ($rmsd < $minRMSD) {
                $minRMSD = $rmsd;
                $min_rotate_angle = [@$rotate_angle];
            }
        }
    };

    my $mingeo;
    while ($keepgoing) {
        if (@$key_center == 2) {
            my $point = $self->get_point($key_center->[0]);
            my $axis;
            if ($key_center->[1] =~ /^\d+$/) {
                $axis = $self->get_bond($key_center->[1], $key_center->[0]);
            }else {
                $axis = $key_center->[1] - $self->get_point($key_center->[0]);
            }

            $self->center_genrotate($point, $axis, deg2rad($angle));
            &__center_genrotate($ref1, $point, $axis, deg2rad($angle));
            $ligand = $self->copy();
            &$rmsd_align([@rotatable_bonds]);

            if ($minRMSD < $last_rmsd ||
                ($wrong_direction < 1)) {
                $keepgoing ++;
                $angle = $direction_changed ? -15 : 15;

                if ($minRMSD < $last_rmsd) {
                    $mingeo = $self->copy();
                    for my $i (0..$#fragments) {
                        my $point = $rotatable_bonds[$i][1];
                        my $axis = $mingeo->get_bond(reverse @{ $rotatable_bonds[$i] });
                        $mingeo->center_genrotate($rotatable_bonds[$i][1], $axis,
                                                  deg2rad($min_rotate_angle->[$i]), $fragments[$i]);
                    }
                    $last_rmsd = $minRMSD;
                }else {
                    $wrong_direction ++;
                }
            }elsif (!$direction_changed) {
                $direction_changed = 1;
                $angle = -1 * $angle * $keepgoing;
                $wrong_direction = 0;
                $keepgoing = 1;
            }else {
                $keepgoing = 0;
                if ($minRMSD > 0.4) {
                    print "Mapping may failed";
                }
            }

            $min_rotate_angle = [(0) x @rotatable_bonds];
            $rotate_ith = -1;

        }else {
            &$rmsd_align($ref1, [@rotatable_bonds], [@fragments]);
            $keepgoing = 0;
        }
    }

    return ([@subs], $mingeo);
}


sub read_geometry {
    my ($self, $file) = @_;

    my ($elements, $flags, $coords, $constraints, $ligand,
        $TM, $key_atoms, $bonds) = AaronTools::FileReader::grab_coords($file);

    $self->{elements} = $elements;
    $self->{flags} = $flags;
    $self->{coords} = $coords;
    $self->{constraints} = $constraints;
    $self->{active_centers} = $key_atoms;

    $self->{RMSD_bonds} = $bonds;

    $self->refresh_connected();
}


#########################################
#some internal function you should never#
#call outside                           #
#########################################

sub __center_genrotate {
    my ($ref, $point, $axis, $degree) = @_;

    my $a = cos($degree/2);
    $axis /= $axis->norm();
    $axis *= sin($degree/2);

    for my $vec (@$ref) {
        if ($vec !~ /^\d+$/) {
            $vec -= $point;

            my $wx = $axis x $vec;
            my $new_vec = $vec + 2*$a*$wx + 2*($axis x $wx);
            map {$vec->[$_] = $new_vec->[$_]} (0..2);
            $vec += $point;
        }
    }
}



package AaronTools::Substrate;
use strict; use warnings;
use Math::Trig;
use Math::Vector::Real;
use Data::Dumper;
our @ISA = qw(AaronTools::Component AaronTools::Geometry);

sub detect_backbone_subs {
    my $self = shift;
    my @sub_atoms;

    for my $target ( sort { $a <=> $b } keys %{ $self->{substituents} } ) {

        my $nearst = $self->{substituents}->{$target}->{end};
        my $groups = $nearst ?
                     $self->get_all_connected($target, $nearst) : '';

        unless ($nearst && $groups) {
            my $min = 999;
            for my $near (@{ $self->{connection}->[$target] }) {
                my $connected = $self->get_all_connected($target, $near);
                if ($#{ $connected } < $min) {
                    $min = $#{ $connected };
                    $nearst = $near;
                    $groups = $connected;
                }
            }
        }

        unless (@{$self->{substituents}->{$target}->{coords}}) {

            my $sub = $self->detect_substituent( target => $target,
                                                    end => $nearst );
            $sub->{sub} = $self->{substituents}->{$target}->{sub};

            $self->{substituents}->{$target} = $sub;
        }

        push (@sub_atoms, @$groups[1..$#{$groups}]);
    }

    my @all_atoms = (0..$#{ $self->{elements} });

    my %all_atoms;
    @all_atoms{ @all_atoms } = ();

    delete @all_atoms{ @sub_atoms };

    my @backbone_atoms = sort { $a <=> $b } keys %all_atoms;

    $self->{backbone} = $self->subgeo([@backbone_atoms]);

    my $old_new = { map { $backbone_atoms[$_] => $_ } (0..$#backbone_atoms) };
    $self->_rearrange_active_con($old_new);
    $self->rebuild();
}


sub copy {
    my $self = shift;

    my $new =  new AaronTools::Geometry();

    bless $new, "AaronTools::Substrate";

    $new->{name} = $self->{name};
    $new->{constraint} = [ map { [ @$_ ] } @{ $self->{constraint} } ];

    if ($self->{substituents}) {
        $new->{substituents} = {};
        for my $target ( sort { $a <=> $b } keys %{ $self->{substituents} } ) {
            $new->{substituents}->{$target} =
                $self->{substituents}->{$target}->copy();
        }
    }

    if ($self->{backbone}) {
        $new->{backbone} = $self->{backbone}->copy();
    }

    if ($self->{constraints}) {
        $new->{constraints} = [ map { [ @$_ ] } @{ $self->{constraints} }];
    }

    if ($self->{backbone}) {
        $new->rebuild();
    }else {
        $new->{elements} = [ @{ $self->{elements} } ];
        $new->{flags} = [ @{ $self->{flags} } ];
        $new->{coords} = [ map { [ @$_ ] } @{ $self->{coords} }];
    }
    $new->refresh_connected();

    return $new;
}




