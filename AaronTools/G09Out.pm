#Contributors: Yanfei Guan and Steven E. Wheeler
#This is a module containing the G09 output class
#This class is initiated from a .log file
use lib $ENV{'QCHASM'};

use AaronTools::Constants qw(PHYSICAL UNIT);
use AaronTools::Atoms qw(:BASIC);

my $elements = ELEMENTS;
my $radii = RADII;

my $h = PHYSICAL->{PLANK};
my $kb = PHYSICAL->{KB};
my $c = PHYSICAL->{SPEED_OF_LIGHT};
my $R = PHYSICAL->{GAS_CONSTANT};
my $P = PHYSICAL->{STANDARD_PRESSURE};
my $hart2kcal = UNIT->{HART_TO_KCAL};
my $amu2kg = UNIT->{AMU_TO_KG};

package AaronTools::G09Out;
use strict; use warnings;
use Math::Trig;
use AaronTools::Geometry;

my $NORM_FINISH = "Normal termination";

my %errors = (
    "NtrErr Called from FileIO" => 'CHK',         #delete
    "Wrong number of Negative eigenvalues" => 'EIGEN',      #opt=noeigen
    "Convergence failure -- run terminated." => 'CONV',    #scf=xqc
    "Erroneous write" => 'QUOTA',               #check quota and alert user; REMOVE error from end of file!
    "Atoms too close" => 'CLASH',               #flag as CLASH
    "The combination of multiplicity" => 'CHARGEMULT',      #die and alert user to check catalyst structure or fix reaction_data!
    "Bend failed for angle" => 'REDUND',                        #Using opt=cartesian
    "Unknown message" => 'UNKNOWN',
);


sub new {
    my $class = shift;
    my %params = @_;

    my $self = {
        file => $params{file},
        read_opt => $params{read_opt},
    };

    bless $self, $class;

    open INFILE, "<$self->{file}" or die "Can't open $self->{file}\n";

    my @coords;
    my @atoms;
    my $finished = 0;
    my $energy = '';
    my $error = '';
    my $error_msg = '';
    my $gradient = {};
    my $frequency;
    my $mass;
    my $T;

    my $hpmodes = 0;
    
    my $Elec_zpve;
    my $enthalpy;
    my $G;
    my $ZPVE;
    my $sigmar;
    my $mult;
    my $charge;
    my @rottemps;
    my @opts;

    while (<INFILE>) {
        my $line = $_;
        /orientation:/  && do  {
            @coords = ();
            @atoms = ();
            $line = <INFILE>;
            $line = <INFILE>;
            $line = <INFILE>;
            $line = <INFILE>;
            $line = <INFILE>;
            do {
                if($line =~ /^\s+\d+\s+(\S+)\s+\S*\s+(\S+)\s+(\S+)\s+(\S+)/) {
                    my $coord = [$2, $3, $4];
                    push(@coords, $coord);
                    push(@atoms, $elements->[$1]);
                }
                $line = <INFILE>;
            } while(!($line =~ /--/));
            if ($self->{read_opt}) {
                push (@opts, [@coords]);
            }
            next;
        };

        /Symbolic Z-matrix:/  && do  {
            @coords = ();
            @atoms = ();
            $line = <INFILE>;
            if ($line =~ /Charge\s=\s*(\d+)\s*Multiplicity\s=\s*(\d+)/) {
                $mult = $2;
                $charge = $1;
            }
            do {
                if($line =~ /^\s*(\S+)\s+([-+]?[0-9]*\.?[0-9]*)\s+([-+]?[0-9]*\.?[0-9]*)\s+([-+]?[0-9]*\.?[0-9]*)/) {
                    my $coord = [$2, $3, $4];
                    push(@coords, $coord);
                    push(@atoms, $1);
                }
                $line = <INFILE>;
            } while(!($line =~ /^\s*$/));
            next;
        };

        /$NORM_FINISH/ && do {
            $finished = 1;
            next;
        };
        #SCF energy
        /SCF Done/o && do {
            my @array = split(/\s+/o, $line);
            $energy = $array[5];
        };

        #Frequencies
        /hpmodes/ && do { $hpmodes = 1; };

        /Harmonic frequencies/ && do {
            my $frequency_file = $_;
            while (<INFILE>) {
                /^\n$/ && do {last;};
                $frequency_file .= $_;
            }

            $frequency = new AaronTools::Frequency( string => $frequency_file,
                                                   hpmodes => $hpmodes,
                                                 num_atoms => scalar @atoms);
        };

        #THERMO
        /Molecular mass:\s+(\S+)/ && do {
            $mass = $1 * $amu2kg;
        };

        /Temperature\s+(\S+)/ && do {
            $T = $1;
        };

        /Rotational constants\s+\(GHZ\):\s+(\S+)\s+(\S+)\s+(\S+)/ && do {
            @rottemps = ($1, $2, $3);
            for my $rot (0..$#rottemps) {
                $rottemps[$rot] *= $h*(10**9)/$kb;
            }
        };

        /Sum of electronic and zero-point Energies=\s+(\S+)/ && do {
            $Elec_zpve = $1;
        };

        /Sum of electronic and thermal Enthalpies=\s+(\S+)/ && do {
            $enthalpy = $1;
        };

        /Sum of electronic and thermal Free Energies=\s+(\S+)/ && do {
            $G = $1;
        };

        /Zero-point correction=\s+(\S+)/ && do {
            $ZPVE = $1;
        };

        /Multiplicity = (\d+)/ && do {
            $mult = $1;
        };

        /Rotational symmetry number\s+(\d+)/ && do {
            $sigmar = $1;
        };
    
        #Gradient
        /Threshold  Converged/ && do {
            while (<INFILE>) {
                /Predicted change in Energy/ && do { last; };
                /(Maximum Force)\s+(\d\.\d+)\s+(\d+\.\d+)\s+(\S+)/ && do {
                    $gradient->{1} = {  item => $1,
                                       thres => $3,
                                       value => $2,
                                   converged => $4 };
                };
                /(RMS     Force)\s+(\d\.\d+)\s+(\d+\.\d+)\s+(\S+)/ && do {
                    $gradient->{2} = {  item => $1,
                                       thres => $3,
                                       value => $2,
                                   converged => $4 };
                };
                /(Maximum Displacement)\s+(\d\.\d+)\s+(\d+\.\d+)\s+(\S+)/ && do {
                    $gradient->{3} = {  item => $1,
                                       thres => $3,
                                       value => $2,
                                   converged => $4 };
                };
                /(RMS     Displacement)\s+(\d\.\d+)\s+(\d+\.\d+)\s+(\S+)/ && do {
                    $gradient->{4} = { item => $1,
                                      thres => $3,
                                      value => $2,
                                  converged => $4 };
                };
            }
        };
       
        (my $temp) = grep { $line =~ /$_/ } keys %errors;
        $error = $errors{$temp} if $temp;
        $error_msg = $line if $temp;
    }

    close(INFILE);

    $self->{geometry} = new AaronTools::Geometry( name => $self->{file},
                                               coords => [@coords],
                                             elements => [@atoms] );
    if (@opts) {
        my $i = 1;
        my @opts_geo;
        for my $coords (@opts) {
            my $new_geo = new AaronTools::Geometry( name => "$self->{file}_$i",
                                                    coords => $coords,
                                                    elements => [@atoms], );
            push (@opts_geo, $new_geo);
        }
        $self->{opts} = [@opts_geo];
    }

    $self->{energy} = $energy;
    $self->{error} = $error;
    $self->{error_msg} = $error_msg;
    $self->{gradient} = $gradient;
    $self->{finished} = $finished;
    $self->{frequency} = $frequency;

    $self->{mass} = $mass;
    $self->{temperature} = $T;

    $self->{rotational_temperature} = [@rottemps];
    $self->{enthalpy} = $enthalpy;
    $self->{free_energy} = $G;

    $self->{multiplicity} = $mult;
    $self->{charge} = $charge;
    $self->{rotational_symmetry_number} = $sigmar;

    if ( (! $self->{finished}) && (! $self->{error} ) ) {
         $self->{error} = $errors{"Unknown message"};
         $self->{error_msg} = "Unknown message";
    }

    return $self;
}


sub geometry {
    my $self = shift;

    return $self->{geometry};
}

sub opts {
    my $self = shift;

    if ($self->{opts}) {
        return @{$self->{opts}};
    }else {
        return ($self->{geometry});
    }
}


sub finished_normal {
    my $self = shift;

    return $self->{finished};
}


sub energy {
    my $self = shift;

    return $self->{energy};
}


sub error {
    my $self = shift;

    return $self->{error};
}


sub error_msg {
    my $self = shift;

    return $self->{error_msg};
}


sub gradient {
    my $self = shift;

    my $gradient = '';

    if ($self->{gradient}) {
        my $ref = $self->{gradient};

        my @gradient = map { "$ref->{$_}->{value}/$ref->{$_}->{converged}" } 
                            keys %{ $ref };

        $gradient = join(", ", @gradient);

        $gradient = $gradient ? "Progress: $gradient" : "Progress: Not found";
    }

    return $gradient;
}


sub frequency {
    my $self = shift;

    return $self->{frequency};
}


sub enthalpy {
    my $self = shift;

    if ($self->{enthalpy}) {
        return $self->{enthalpy};
    }else {
        print "No enthalpy found\n";
    }
}


sub free_energy {
    my $self = shift;

    if ($self->{free_energy}) {
        return $self->{free_energy};
    }else {
        print "No free energy found\n";
    }
}


sub Grimme_G {
    my $self = shift;

    my %params = @_;

    my $v0 = 100;           #cutoff for quasi rrho

    if (!$self->frequency()) {
        print "Cannot calculate the Grimme free energy without frequncy information\n";
    }
    my $rottemps = $self->{rotational_temperature};

    my $Bav = ($h**2/(24*pi**2*$kb))*(1/$rottemps->[0] + 1/$rottemps->[1] + 1/$rottemps->[2]);

    my $T = $params{temperature} || $self->{temperature};
    my $mass = $self->{mass};
    my $sigmar = $self->{rotational_symmetry_number};
    my $mult = $self->{multiplicity};

    my @freqs = $self->frequency()->positive_frequencies();
    my @vibtemps = map { $_*$c*$h/$kb } @freqs;

    #Translational
    my $qt = (2*pi*$mass*$kb*$T/($h*$h))**(3/2)*$kb*$T/$P;
    my $St = $R*(log($qt) + 5/2);
    my $Et = 3*$R*$T/2;

    #Electronic
    my $Se = $R*(log($mult));

    #Rotational
    my $qr = (sqrt(pi)/$sigmar)*($T**(3/2)/sqrt($rottemps->[0]*$rottemps->[1]*$rottemps->[2]));
    my $Sr = $R*(log($qr) + 3/2);
    my $Er = 3*$R*$T/2;

    #Vibirational
    my ($Ev, $Sv_quasiRRHO) = (0, 0, 0);

    for my $i (0..$#freqs) {
        my $Sv_temp = $vibtemps[$i]/($T*(exp($vibtemps[$i]/$T)-1)) - 
                      log(1-exp(-$vibtemps[$i]/$T));
        $Ev += $vibtemps[$i]*(1/2 + 1/(exp($vibtemps[$i]/$T) - 1));

        #quassi_RRHO
        my $mu = $h/(8*pi**2*$freqs[$i]*$c);
        my $mu_prime = $mu*$Bav/($mu + $Bav);
        my $Sr_eff = 1/2 + log(sqrt(8*pi**3*$mu_prime*$kb*$T/$h**2));

        my $weight = 1/(1+($v0/$freqs[$i])**4);

        $Sv_quasiRRHO += $weight*$Sv_temp + (1-$weight)*$Sr_eff;
    }

    $Ev *= $R; $Sv_quasiRRHO *= $R;

    my @a = ($St, $Sr, $Sv_quasiRRHO, $Se);

    my $Hcorr = $Et + $Er + $Ev + $R*$T; 
    my $Stot_quasiRRHO = $St + $Sr + $Sv_quasiRRHO + $Se;
    my $Gcorr_quasiRRHO = $Hcorr - $T*$Stot_quasiRRHO;
    $Gcorr_quasiRRHO /= $hart2kcal * 1000;

    my $Grimme_G = $self->{energy} + $Gcorr_quasiRRHO;

    return $Grimme_G;
}

sub bond_change {
    my $self = shift;
    my ($bond) = @_;

    my $d_ref = $self->{opts}->[0]->distance( atom1=>$bond->[0], atom2=>$bond->[1] );

    my $n=$#{$self->{opts}};
    for (my $i=$#{$self->{opts}}; $i>=0; $i--) {
        my $d = $d_ref - 
            $self->{opts}->[$i]->distance(atom1=>$bond->[0], atom2=>$bond->[1]);
        if (abs($d) < 0.25) {
            $n = $i;
            last;
        }
    }

    return $n;
}



package AaronTools::Frequency;
use strict; use warnings;
use AaronTools::Atoms qw(:BASIC);


sub new {
    my $class = shift;
    my %params = @_;

    my ($string, $hpmodes, $num_atoms) = ( $params{string}, $params{hpmodes},
                                           $params{num_atoms} );
    my @input = split(/\n/, $string);

    my @num_head = grep { $input[$_] =~ /Harmonic frequencies/ } (0..$#input);

    if ($hpmodes && (@num_head == 2)) {
        @input = @input[0..$num_head[1] - 1];
    }elsif ($hpmodes) {
        print "The .log file is damaged, cannot get frequency\n";
        return;
    }

    my @frequencies;
    my @intensities;
    my @vectors;
    my $num_coords;
    my $vectors_block;

    while(@input) {
        my $line = shift @input;

        if ($line =~ /Frequencies/) {
            
            $line =~ s/^\s+//;

            my @patterns = split(/\s+/, $line);
            shift @patterns; shift @patterns;

            if ($patterns[0] !~ /^[-+]?[0-9]*\.?[0-9]*$/) {
                my $string = shift @patterns;
                my @patterns_temp = unpack("(A10)*", $string);

                unshift(@patterns, @patterns_temp);
            }

            push (@frequencies, @patterns);
            shift @input, shift @input;
            next;
        }

        if ($line =~ /IR Inten/) {

            $line =~ s/^\s+//;

            my @patterns = split(/\s+/, $line);

            shift @patterns; shift @patterns; shift @patterns;

            push (@intensities, @patterns);
            shift @input;

            $num_coords = 0;
            $vectors_block = [];
            next;
        }
        
        #vectors
        if ($hpmodes) {
            if ($line =~ 
                /\s+(\d+)\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {

                my @moves = ($3, $4, $5, $6, $7);
                my $atom = $2;
                for my $i (0..4) {
                    $vectors_block->[$i]->[$atom-1]->[$1-1] = $moves[$i];
                }
                $num_coords ++;

                if ($num_coords == $num_atoms * 3) {
                    #deep copy of vectors
                    push(@vectors, map { [map { [@$_] } @$_] } @$vectors_block);
                }
                next;
            }
        }else{
            if ($line =~
                /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {

                my $moves = [ [ $2, $3, $4 ],
                              [ $5, $6, $7 ],
                              [ $8, $9, $10] ];

                my $atom = $1;

                for my $i (0..2) {
                    for my $xyz (0..2) {
                        $vectors_block->[$i]->[$atom-1]->[$xyz] = $moves->[$i]->[$xyz];
                    }
                }
                
                if ($atom == $num_atoms) {
                    push (@vectors, map { [map { [@$_] } @$_] } @$vectors_block);
                }
                next;
            }
        }
    }

    my $self = { 
        frequencies => [@frequencies],
        intensities => [@intensities],
            vectors => [@vectors],
    };

    bless $self, $class;

    return $self;
}
        

sub frequencies {
    my $self = shift;

    return @{ $self->{frequencies} };
}


sub imaginary_frequencies {
    my $self = shift;

    my @imaginary = grep { $_ < 0 } $self->frequencies();

    return @imaginary;
}


sub positive_frequencies {
    my $self = shift;

    my @p_frequencies = grep { $_ > 0 } $self->frequencies();

    return @p_frequencies;
}


sub lowest_frequency {
    my $self = shift;

    return $self->{frequencies}->[0];
}


sub intensity {
    my $self = shift;

    my ($frequency) = @_;

    my @frequencies = $self->frequencies();

    my ($num) = grep { $frequencies[$_] == $frequency } (0..$#frequencies);

    return $self->{intensities}->[$num];
}


sub vector {
    my $self = shift;

    my ($frequency) = @_;

    my @frequencies = $self->frequencies();

    my ($num) = grep { $frequencies[$_] == $frequency } (0..$#frequencies);

    return $self->{vectors}->[$num];
}


sub is_TS {
    my $self = shift;

    if ($self->imaginary_frequencies() == 1) {
        return 1;
    }else { return 0;}
}























    






