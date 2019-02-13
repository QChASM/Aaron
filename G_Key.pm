use lib $ENV{'QCHASM'};
use AaronTools::Constants qw(NAMES);

my $TS_lib = NAMES->{TS_LIB};

package Aaron::G_Key;
use strict; use warnings;
use Data::Dumper;

my $HOME = $ENV{'HOME'};
my $QCHASM = $ENV{'QCHASM'};
$QCHASM =~ s|/\z||;	#Strip trailing / from $QCHASM if it exists

sub new {
    my $class = shift;
    my %params = @_;
    my $self = {
        solvent => $params{solvent},            #solvent (gas, dichloromethane, etc)
        temperature => $params{colvent},        #temperature
        denfit => $params{denfit},              #Boolean for whether to use denfit
        charge => $params{charge},              #charge
        mult => $params{mult},                  #multiplicity
        solvent_model => $params{solvent_model},#Solvent model (pcm, smd, etc)
        n_procs => $params{n_procs},            #number of cores
        wall => $params{wall},                  #walltime
        short_procs => $params{short_procs},    #number of cores for short jobs
        short_wall => $params{short_wall},      #walltime for short jobs
        node=> $params{node},                   #Node type?
        grid=> $params{grid},                   #integration grid keyword
        gen => $params{gen},                    #gen basis
        custom => $params{custom},              #Custom default from .aaronrc
        emp_dispersion => $params{emp_dispersion},          #Empirical dispersion
        con_thres => $params{con_thres},        #connectivity threshold
    };

    bless $self, $class;

    $self->{level} = new Aaron::Theory_level();
    $self->{low_level} = new Aaron::Theory_level();
    $self->{high_level} = new Aaron::Theory_level();

    return $self;
}


sub read_key_from_com {
    my $self = shift;

    my ($com, $level) = @_;

    $level //= 'level';

    open COM, "<$com" or die "Cannot open $com:$!\n";

    my $command_line = <COM>;

    if ($command_line =~ /^\%/) {
        $command_line = <COM>;
    }

    my $temperature;
    my $solvent;
    my $solvent_model;
    ($self->{$level}->{final_method}) = $command_line =~ /^\#(\S+)(\s+)/;
    ($temperature) = $command_line =~ /temperature=(\d+\.?\d+)/;
    ($solvent) = $command_line =~ /solvent=(\S+\))/;
    ($solvent_model) = $command_line =~ /scrf=\((\S+),/;
    ($self->{emp_dispersion}) = $command_line =~ /EmpiricalDispersion=(\S+)/;
    ($self->{grid}) = $command_line =~ /int=\(grid=(\S+)/;

    $temperature && do {$self->{temperature} = $temperature};
    $solvent && do {$self->{solvent} = $solvent};
    $solvent_model && do {$self->{solvent_model} = $solvent_model};

    my $chargemult;

    for (1..4) {
        $chargemult = <COM>;
    }

    ($self->{charge}, $self->{mult}) = $chargemult =~ /^(\d)\s(\d)/;


    my $hit_coords;
    while (<COM>) {
        /^\s?\S+\s+[+-]?\d+\.?\d+\s+[+-]?\d+\.?\d+\s+[+-]?\d+\.?\d+/ && do {
            $hit_coords = 1;
        };
        /^$/ && do {last if $hit_coords};
    }

    my $footer = '';
    my $hit_constraint;
    while (<COM>) {
        /^B\s\d+\s\d+\sF/ && do { 
            while (<COM>) {
                /^$/ && do {last};
            }
            next;
        };
        $footer .= $_;
    }

    $self->{$level}->{footer} = $footer;

    close COM;
}


sub _read_key_from_input {
    my $self = shift;

    my ($input) = @_;

    my $theory = $self->{custom};

    open INPUT, "<$input" or die "Can't open $input:$!\n";

    my $hit = $theory ? 0 : 1;

    my $read_low_level = $self->{low_level}->method() ? 0 : 1;
    my $read_level = $self->{level}->method() ? 0 : 1;
    my $read_high_level = $self->{high_level}->method() ? 0 : 1;

    while (<INPUT>) {
        if ($theory) {
            /^$theory$/ && do {$hit = 1; next};
        }

        /^$/ && do {last if $hit};
        /[Cc]ustom=(\S+)/ && do {$self->{custom} = $1; next;};
        /^\s*[gG]en=(\S+)/ && do {$self->{gen} = $1 unless $self->{gen}; next;};

        if ($hit) {
            #system
            /^\s*[nN]_procs=(\d+)/ && do{
                $self->{n_procs} = $1 unless $self->{n_procs}; next;
            };
            
            /^\s*[wW]all=(\S+)/ && do {
                $self->{wall} = $1 unless $self->{wall}; next;};
            /^\s*[sS]hort_procs=(\d+)/ && do {
                $self->{short_procs} = $1 unless $self->{short_procs}; next;
            };

            /^\s*[sS]hort_wall=(\S+)/ && do {
                $self->{short_wall} = $1 unless $self->{short_wall}; next;
            };

            /^\s*[nN]ode=(\s)/ && do {
                $self->{node} = $1 unless $self->{node}; next;
            };
            #G09
            /grid=(\S+)/ && do {$self->{grid} = $1 unless $self->{grid}; next;};
            /\b[sS]olvent=(\S+)/ && do {$self->{solvent} = $1 unless $self->{solvent}; next;};
            /\b[sS]olvent_model=(\S+)/ && do {$self->{solvent_model} = $1 unless $self->{solvent_model}; next;};
            /\b[hH]igh_solvent=(\S+)/ && do {$self->{high_solvent} = $1 unless $self->{high_solvent}; next;};
            #input file line for using a solvent for the high-level single point is high_solvent=X
            /\b[hH]igh_solvent_model=(\S+)/ && do {$self->{high_solvent_model} = $1 unless $self->{high_solvent_model}; next;};
            #solvent model for high-level single point
            /\s*[tT]emperature=(\S+)/ && do {
                $self->{temperature} = $1 unless $self->{temperature}; next;
            };

            /^\s*[eE]mp_dispersion=(\S+)/ && do {
                $self->{emp_dispersion} = $1 unless $self->{emp_dispersion}; next;
            };

            if ($read_low_level) {
                /^\s*[lL]ow_method=(\S+)/ && do {$self->{low_level}->read_method($1); next;};
                /^\s*[lL]ow_basis=(.+)/ && do {$self->{low_level}->read_basis($1); next;};
                /^\s*[lL]ow_ecp=(\S+)/ && do {$self->{low_level}->read_ecp($1); next;};
            }

            if ($read_level) {
                /^\s*[mM]ethod=(\S+)/ && do {$self->{level}->read_method($1); next;};
                /^\s*[bB]asis=(.+)/ && do {$self->{level}->read_basis($1); next;};
                /^\s*[eE]cp=(.+)/ && do {$self->{level}->read_ecp($1); next;};
            }

            if ($read_high_level &&
                (!$theory || ($theory ne 'Default'))) {
                /^\s*[hH]igh_method=(\S+)/ && do {$self->{high_level}->read_method($1); next;};
                /^\s*[hH]igh_ecp=(.+)/ && do {$self->{high_level}->read_ecp($1); next;};
                /^\s*[hH]igh_basis=(.+)/ && do {$self->{high_level}->read_basis($1); next;};
            }

            /\s*[dD]enfit=(\S+)/ && do {$self->{denfit} = $1 unless defined $self->{denfit}; next;};
            /\s*[cC]harge=(\S+)/ && do {$self->{charge} = $1 unless defined $self->{charge}; next;};
            /\s*[mM]ult=(\S+)/ && do {$self->{mult} = $1 unless defined $self->{mult}; next;}; 
            /\s*[cC]on_thres=(\S+)/ && do {$self->{con_thres} = $1 unless defined $self->{con_thres}; next;};
        }
        
    }

    close INPUT;
    return $hit;
}


sub read_key_from_input {
    my $self = shift;

    my ($input) = @_;
 
    $self->_read_key_from_input($input) if $input;

    $self->{custom} //= 'Default';

#This should read global defaults from $QCHASM/Aaron/.aaronrc and then overwrite any of these with personal defaults from $HOME/.aaronrc
#This way people can just set their own defaults for the keywords that are different
    my $hit;
    if (-e "$HOME/.aaronrc") {
        $self->_read_key_from_input("$HOME/.aaronrc");
    }

    if (-e "$QCHASM/Aaron/.aaronrc") {
        $self->_read_key_from_input("$QCHASM/Aaron/.aaronrc");
    }

    for my $level ( $self->{low_level}, $self->{level}, $self->{high_level} ) {
        $level->check_gen($self->{gen});
    }

    $self->{solvent} //= 'gas';
    #solvent defaults to gas if none is given
    $self->{high_solvent} //= $self->{solvent};
    #high_solvent defaults to solvent if none is given
    $self->{con_thres} //= 0.5;

    $self->{solvent_model} //= $self->{solvent} =~ /^[Gg]as$/ ? '' : 'pcm'; 
    #solvent model defaults to PCM if a solvent is given
    $self->{high_solvent_model} //= $self->{high_solvent} =~ /^[Gg]as$/ ? '' : 'pcm'; 
}
        

package Aaron::Workflow_Key;
use strict; use warnings;
use Data::Dumper;

sub new {
    my $class =shift;
    my %params = @_;
    my $self = {
        reaction_type => $params{reaction_type},
        template => $params{template},
        step => $params{step},
        input_conformers_only => $params{input_conformers_only},
        full_conformers => $params{full_conformers},
        selectivity => $params{selectivity},
        multistep => $params{multistep},
        debug => $params{debug},
        nosub => $params{nosub},
        record => $params{record},
        short => $params{short},
        no_quota => $params{no_quota},
   };

   bless $self, $class;

   $self->{reaction_type} //= '';
   $self->{template} //= '';
   $self->{step} //= [];
   $self->{input_conformers_only} //= 0;
   $self->{full_conformers} //= 0;
   $self->{selectivity} //= ['R', 'S'];
   $self->{multistep} //= 0;
   $self->{debug} //= 0;
   $self->{nosub} //= 0;
   $self->{record} //= 0;
   $self->{short} //= 0;
   $self->{no_quota} //= 0;
   
   return $self;
}

sub read_input {
    my $self = shift;

    my ($input) = @_;

    open IN, "<$input" or die "Cannot open $input:$!\n";

    while(<IN>) {
        /[rR]eaction_type=(\S+)/ && do {$self->{reaction_type} = $1; next;};
        /\b[tT]emplate=(\S+)/ && do {$self->{template} = $1; next;};
        /[sS]tep=(\S+)/ && do {$self->{step} = [split(/;/, $1)]; next;};
        /[Ii]nput_conformers_only=(\d+)/ && do {$self->{input_conformers_only} = $1; next;};
        /[Ff]ull_conformers=(\d+)/ && do {$self->{full_conformers} = $1; next;};
        /[Ss]electivity=(\S+)/ && do {$self->{selectivity} = [split(/;/, $1)];next;};
        /[Mm]ultistep=(\d+)/ && do {$self->{multistep} = $1; next;};
        #non-default job submission template
        /[sS]ubmission_template=(\S+)/ && do {$self->{submission_template} = $1; next; };
    }

    close IN;
    $self->examine();
}

sub examine {
    my $self = shift;

    #chomp each value
    for my $key (keys %{ $self }) {
        chomp($self->{$key});
    }
    
    #examine arguments;
    if ($#{$self->{selectivity}} == 0 &&
        ($self->{selectivity}->[0] eq 'NONE')) {
        $self->{selectivity} = [];
    }
    
    unless ($self->{template}) {
        print "A template must be given explicitly in the <jobname>.in file " .
              "by template=xxxxx, \n" .
              "If this catalyst contains different steps in a reaction, " .
              "you should also specify the step too. e.g. template=catalyst/TS1. \n" .
              "Exiting without running any calculations\n";
        exit(1);
    }

    my $TS_path = (-d "$HOME/Aaron_libs/$TS_lib/$self->{reaction_type}/$self->{template}") ?
                    "$HOME/Aaron_libs/$TS_lib/$self->{reaction_type}/" : 
                    "$QCHASM/Aaron/$TS_lib/$self->{reaction_type}/";

    unless (-d $TS_path) {
        print "Can't find $self->{template} in either user defined TS library: ".
              "$HOME/Aaron_libs/$TS_lib/ or the built_in library: $QCHASM/Aaron/$TS_lib/$self->{template}\n";
        exit(1);
    }

    if (@{ $self->{step} } > 1) {
        $self->{multistep} = 1;
    }

    $self->{TS_path} = $TS_path;
}


package Aaron::Theory_level;
use strict; use warnings;
use lib $ENV{'QCHASM'};
use Data::Dumper;

use AaronTools::Atoms qw(:BASIC);

my $masses = MASS;
my $tmetal = TMETAL;

sub new {
    my $class = shift;

    my %params = @_;

    my $self = { 
        basis => {},
        gen_basis => [],
        method => $params{method},
    };

    bless $self, $class;

    return $self;
}


sub read_method {
    my $self = shift;

    my ($method) = @_;

    $self->{method} = $method;
}


sub read_basis {
    my $self = shift;

    my ($line) = @_;

    my @entry = split(/\s+/, $line);

    my @atoms;
    for my $entry (@entry) {
        if (exists $masses->{$entry}) {
            push (@atoms, $entry);
        }elsif ($entry =~ /^tm$/i) {
            for my $element (keys %{ $tmetal }) {
                push (@atoms, $element);
            }
        }else {
            if ($entry =~ /^gen/) {
                push (@{$self->{gen_basis}}, $entry);
            }elsif (@atoms) {
                @{$self->{basis}}{@atoms} = ($entry) x (scalar @atoms);
            }else {
                for my $element (keys %{ $masses }) {
                    unless (exists $self->{basis}->{$element}) {
                        $self->{basis}->{$element} = $entry;
                    }
                }
            }
        }
    }
}


sub read_ecp {
    my $self = shift;

    my ($line) = @_;

    $self->{ecp} = {};

    my @entry = split(/\s+/, $line);

    my @atoms;
    for my $entry (@entry) {
        if (exists $masses->{$entry}) {
            push (@atoms, $entry);
        }elsif ($entry =~ /^tm$/i) {
            for my $element (keys %{ $tmetal }) {
                push (@atoms, $element);
            }
        }else {
            if (@atoms) {
                @{$self->{ecp}}{@atoms} = ($entry) x (scalar @atoms);
            }else {
                print "You must indicate the atom type using this ecp: $entry\n";
                exit 0;
            }
        }
    }
}


sub check_gen {
    my $self = shift;

    my ($gen) = @_;

    $gen //= '';

    my @gen_basis = @{ $self->{gen_basis} };
    my @gen_files = grep { $_ ne 'gen' } @gen_basis;

    my @correct_gen_basis;
    my $typein;
    if (@gen_files) {
        if ($gen) {
            my @basis;

            opendir(GEN, $gen) or do {
                my $msg = "Cannot open directory $gen to read basis set files\n".
                          "please stop at this time and fix your gen library\n".
                          "Or type in your basis following the guide.\n";
                warn($msg);
                $typein = 1;
            };
            unless ($typein) {
                $self->{gen} = $gen;
                $self->{gen} =~ s/\/$//;
                $self->{gen} .= '/';

                @basis = readdir(GEN);
                closedir(GEN);
                for my $file (@gen_files) {
                    my $basis_file = $file;
                    $basis_file =~ s/^gen\///;

                    unless (grep { $_ eq $basis_file } @basis) {
                        my $msg = "Cannot find basis set file $basis_file in $gen\n".
                                  "please stop at this time and fix your gen library\n".
                                  "Or type in your basis following the guide.\n";
                        warn($msg);
                        $typein = 1;
                    }else {
                        push (@correct_gen_basis, $file);
                    }
                    
                    for (keys %{ $self->{basis} }) {
                        delete $self->{basis}->{$_} if ($self->{basis}->{$_} eq $file);
                    }
                }
            }
        }else {
            my $msg = "No path to the basis set files was found.\n".
                      "Please stop at this time and write your path to the basis library in xxxx\n".
                      "Or type in your basis following the guide.\n";
            warn($msg);
            $typein = 1;
        }
    }

    $self->{gen_basis} = [@correct_gen_basis];

    $typein = 1 if grep { $_ eq 'gen' } @gen_basis;

    if ($typein) {
        print "Type in your basis below. For example:\n".
              "-H     0\n".
              "S   3   1.00\n".
              "     34.0613410              0.60251978E-02\n".
              "      5.1235746              0.45021094E-01\n".
              "      1.1646626              0.20189726\n".
              "S   1   1.00\n".
              "      0.32723041             1.0000000\n".
              "S   1   1.00\n".
              "      0.10307241             1.0000000\n".
              "P   1   1.00\n".
              "      0.8000000              1.0000000\n".
              "****\n".
              "-C     0\n".
              "......\n\n".
              "And then press ctrl-D to finishi entering\n";

        my @basis = <STDIN>;
        chomp @basis;

        my $basis = join ("\n", @basis);

        $self->{typein_basis} = $basis;
    }
}


sub method {
    my $self = shift;

    my $method;

    if ($self->{final_method}) {
        return $self->{final_method};
    }

    if ($self->{ecp}) {
        $method = $self->{method} . "/genecp";
    }elsif (@{$self->{gen_basis}}) {
        $method = $self->{method}. "/gen";
    }elsif (@{$self->unique_basis()} > 1) {
        $method = $self->{method}. "/gen";
    }elsif ($self->{typein_basis}) {
        $method = $self->{method}. "/gen";
    }elsif (@{$self->unique_basis()} == 1) {
        $method = $self->{method} . "/" . $self->unique_basis()->[0];
    }else {
        $method = $self->{method};
    }

    return $method;
}

    
sub unique_basis {
    my $self = shift;

    my @basis = values %{ $self->{basis} };

    my @unique;

    for my $basis (@basis) {
        unless (grep {$_ eq $basis} @unique) {push (@unique, $basis);}
    }

    return [@unique];
}


sub unique_ecp {
    my $self = shift;

    my @basis = values %{ $self->{ecp} };

    my @unique;

    for my $basis (@basis) {
        unless (grep {$_ eq $basis} @unique) {push (@unique, $basis);}
    }

    return [@unique];
}

sub footer {
    my $self = shift;
    my ($geometry) = @_;

    my $return = '';

    if ($self->{footer}) {
        return $self->{footer};
    }

    if (keys %{ $self->{basis} }) {
        if (@{$self->unique_basis()} > 1 ||
            $self->{ecp} ||
            @{$self->{gen_basis}} ||
            $self->{typein_basis}) {

            my @elements = @{ $geometry->{elements} };

            my @unique_basis = @{$self->unique_basis()};

            for my $basis (@unique_basis) {
                my @atoms = grep { $self->{basis}->{$_} eq $basis } keys %{ $self->{basis} };

                my @exit_atoms;

                for my $atom (@atoms) {
                    if (grep { $_ eq $atom } @elements) {
                        push (@exit_atoms, $atom);
                    }
                }

                next unless @exit_atoms;

                $return .= sprintf "%s " x @exit_atoms, @exit_atoms;
                $return .= "0\n";
                $return .= "$basis\n";
                $return .= '*' x 4 . "\n";
            }
        }
    }

    if (@{ $self->{gen_basis} }) {
        for my $gen_basis (@{ $self->{gen_basis} }) {
            $return .= "\@" . $self->{gen} . $gen_basis . "/N\n";
        }
    }

    if ($self->{typein_basis}) {
        $return .= "$self->{typein_basis}\n";
    }

    if ($self->{ecp}) {
        $return .= "\n";

        my @elements = @{ $geometry->{elements} };

        my @unique_ecp = @{$self->unique_ecp()};

        for my $ecp (@unique_ecp) {
            my @atoms = grep { $self->{ecp}->{$_} eq $ecp } keys %{ $self->{ecp} };

            my @exit_atoms;

            for my $atom (@atoms) {
                if (grep { $_ eq $atom } @elements) {
                    push (@exit_atoms, $atom);
                }
            }

            next unless @exit_atoms;

            $return .= sprintf "%s " x @exit_atoms, @exit_atoms;
            $return .= "0\n";
            $return .= "$ecp\n";
        }
    }

    return $return;
}

sub footer_log {
    my $self = shift;

    my $return = '';
    my @unique_basis = $self->unique_basis();

    for my $basis (@unique_basis) {
        $return .= "$basis\n";
    }

    $return .= "$self->{gen_basis}\n" if @{$self->{gen_basis}};
    $return .= "$self->{typein_basis}\n" if $self->{typein_basis};

    return $return;
}









1
