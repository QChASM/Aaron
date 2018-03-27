package AaronInit;

use strict; use warnings;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

my $HOME = $ENV{'HOME'};
my $AARON = $ENV{'AARON'};

use Constants qw(:INFORMATION :THEORY :PHYSICAL :SYSTEM :JOB_FILE :OTHER_USEFUL);
use Pod::Usage;

my @authors = @{ INFO->{AUTHORS} };
my $version = INFO->{VERSION};
my $lastupdate = INFO->{LASTUPDATE};

my $module = {};
#Check necessary modules are installed
my $helpMsg = "\nAARON: An Automated Reaction Optimizer for new catalyst\n".
              "- Computational toolkit to aid in optimizing transition states with the Gaussian09 quantum chemistry software.\n\n".
              "Authors: @authors.\n\nLast Update: $lastupdate\n\n".
              "AARON automates the optimization of transition states for a wide variety of metal-free asymmetric reactions.\n".
              "Based on a library of TS structures previously computed using model catalysts, \n".
              "AARON replaces the model catalyst with a user-supplied catalyst and then performs a prescribed series of \n".
              "constrained and uncontrained optimizations to arrive at final predicted structures and energies for all transition states.\n";

#arguments for AARON taking from command line
our %arg_parser = ( sleeptime => SLEEP_TIME );
&read_args();

&check_modules();

use Cwd qw(getcwd);
use Getopt::Long;
use Exporter qw(import);
use AaronTools::Atoms qw(:BASIC :LJ);
use Data::Dumper;

#default values for some argument

our @EXPORT = qw(%arg_in %arg_parser $ligs_subs $parent $jobname
                 $template_job $system init_main grab_cata_coords write_status);

my $TS_lib = NAMES->{TS_LIB};
my $LIG_OLD = NAMES->{LIG_OLD};
my $LIG_NONE = NAMES->{LIG_NONE};
our $jobname;
my $masses = MASS;
our $parent = getcwd();
my $input_file;

#content of template job file
our $template_job = {};

#ligand and substituent information
our $ligs_subs = {};

#all arguments for AARON
our %arg_in = (
    TS_path => '',
    MaxRTS => 1,
    MaxSTS => 1,
    solvent => 'gas',
    temperature => ROOM_TEMPERATURE,
    reaction_type => '',
    gen => '',
    low_level => new Theory_level( method => 'PM6' ),
    level => new Theory_level( method => 'B3LYP' ),
    high_level => new Theory_level(),
    denfit => 1,
    pcm => PCM,
    catalyst => '',
    charge => 0,
    mult => 1,
    template => '',
    emp_dispersion => '',
    input_conformers_only => 0,
    full_conformers => 0,
    selectivity => ['R', 'S'],
    no_ee => 0,
);

#THe default system is Ada if not specified
#$system is a hash reference
our $system = ADA;

#read command line arguments
sub check_modules {
    print "Checking required modules...\n";

    eval {
        require Math::Vector::Real;
        $module->{Real}->{install} = 1;
        1;
    } or do {$module->{Real}->{install} = 0;};
    $module->{Real}->{mod} = 'Math::Vector::Real';

    eval {
        require Math::MatrixReal;
        $module->{Matrix}->{install} = 1;
        1;
    } or do {$module->{Matrix}->{install} = 0;};
    $module->{Matrix}->{mod} = 'Math::MatrixReal';

    eval {
        require Data::Dumper;
        $module->{Dumper}->{install} = 1;
        1;
    } or do {$module->{Dumper}->{install} = 0;};
    $module->{Dumper}->{mod} = 'Data::Dumper';

    my $msg = '';
    my $quit = 0;
    for my $mod (sort keys %$module) {
        if (! $module->{$mod}->{install}) {
            $msg .= "No $module->{$mod}->{mod} found in the perl library. "  .
                    "This is fatal to Aaron, Aaron will quit at this time. " .
                    "Please install that module.\n";
            $quit = 1;
        }
    }
    print $msg;

    exit (1) if $quit;
    
    print "Necessary modules found, Aaron will start in a second\n"; 
}


        

sub read_args{
    GetOptions(
        'debug' => \$arg_parser{debug},
        'nosub' => \$arg_parser{nosub},
        'help|h' => \$arg_parser{help},
        'restart' => \$arg_parser{restart},
        'record' => \$arg_parser{record},
        'short' => \$arg_parser{short},
        'multistep' => \$arg_parser{multistep},
        'absthermo' => \$arg_parser{absthermo},
        'sleep=s' => \$arg_parser{sleeptime},
    ) or pod2usage (
        -input => "$AARON/pod_ref",
        -exitval => 1,
        -verbose => 1 );

    pod2usage(
        -msg => $helpMsg,
        -input => "$AARON/pod_ref",
        -exitval => 1,
        -verbose => 1) if $arg_parser{help};

    ($input_file) = grep { $_ =~ /\.in$/ } @ARGV;

    $input_file or pod2usage (
        -msg => "A input file must be provided\n",
        -input => "$AARON/pod_ref",
        -exitval => 1,
        -verbose => 0 );
}


#read arguments from input file
sub read_params {
    open my $in_h, "< $input_file" or die "Can't open $input_file:$!\n";

    ($jobname) = $input_file =~ /(\S+)\.in/; 

    my $lig = {};
    my $sub = {};

    while(<$in_h>) {
        /[sS]olvent=(\S+)/ && do {$arg_in{solvent} = $1; next;};
        /[pP]cm=(\S+)/ && do {$arg_in{pcm} = $1; next;};
        /[tT]emperature=(\S+)/ && do {$arg_in{temperature} = $1; next;};
        /[rR]eaction_type=(\S+)/ && do {$arg_in{reaction_type} = $1; next;};
        /^[eE]mp_dispersion=(\S+)/ && do {$arg_in{emp_dispersion} = $1; next;};
        /[gG]en=(\S+)/ && do {$arg_in{gen} = $1; next;};
        /[lL]ow_method=(\S+)/ && do {$arg_in{low_level}->read_method($1); next;};
        /^[mM]ethod=(\S+)/ && do {$arg_in{level}->read_method($1); next;};
        /[hH]igh_method=(\S+)/ && do {$arg_in{high_level}->read_method($1); next;};
        /^[lL]ow_ecp=(\S+)/ && do {$arg_in{low_level}->read_ecp($1); next;};
        /^[eE]cp=(.+)/ && do {$arg_in{level}->read_ecp($1); next;};
        /[hH]igh_ecp=(.+)/ && do {$arg_in{high_level}->read_ecp($1); next;};
        /^[lL]ow_basis=(.+)/ && do {$arg_in{low_level}->read_basis($1); next;};
        /^[bB]asis=(.+)/ && do {$arg_in{level}->read_basis($1); next;};
        /^[hH]igh_basis=(.+)/ && do {$arg_in{high_level}->read_basis($1); next;};
        /[dD]enfit=(\S+)/ && do {$arg_in{denfit} = $1; next;};
        /[cC]harge=(\S+)/ && do {$arg_in{charge} = $1; next;};
        /[mM]ult=(\S+)/ && do {$arg_in{mult} = $1; next;}; 
        /[tT]emplate=(\S+)/ && do {$arg_in{template} = $1; next;};
        /[Ii]nput_conformers_only=(\d+)/ && do {$arg_in{input_conformers_only} = $1; next;};
        /[Ff]ull_conformers=(\d+)/ && do {$arg_in{full_conformers} = $1; next;};
        /[Ss]electivity=(\S+)/ && do {$arg_in{selectivity} = [split(/;/, $1)];next;};
        /[Nn]o_ee=(\d+)/ && do {$arg_in{no_ee} = $1; next;};

        #System information
        /^\&[sS]ystem/ && do {
            while(<$in_h>) {
                /^\&$/ && do {last;};
                /n_procs=(\d+)/ && do {$system->{N_PROCS}=$1; next;};
                /^wall=(\d+)/ && do {$system->{WALL}=$1; next;};
                /short_procs=(\d+)/ && do {$system->{SHORT_PROCS}=$1; next;};
                /^short_wall=(\d+)/ && do {$system->{SHORT_WALL}=$1; next;};
                /node_type=(\S+)/ && do {$system->{NODE_TYPE}=$1; next;};
            }
        };
        #new ligand information
        /^\&[Ll]igands/ && do {
            while(<$in_h>) {
                /^\&$/ && do {last;};
                /^\s*(\S+)\:\s*(.*)/ && do {
                    my $lig_ali = $1;
                    my $lig_info = $2;
                    my $lig_new;
                    my $lig_sub;
                    if ($lig_info =~ /^none$/i) {
                        $lig_new = $LIG_NONE;
                    }elsif ($lig_info =~ /^[Ll]igand=(\S+)(\s+(.*))?/) {
                        $lig_new = $1;
                        $lig_sub = $3;
                    }else{ 
                        $lig_new = $LIG_OLD;
                        $lig_sub = $lig_info;
                    }

                    $lig_sub //= '';

                    $lig->{$lig_ali} = {};

                    $lig->{$lig_ali}->{ligand} = $lig_new;
                    my %substituents = split(/[=\s]/, $lig_sub);
                    $lig->{$lig_ali}->{substituents} = { map { $_ - 1 => $substituents{$_} }
                                                            keys %substituents };
                }
            }
        };
        #new Substrate information
        /^\&[Ss]ubstrates/ && do {
            while(<$in_h>) {
                /^\&$/ && do {last;};
                /^(Sub\d+)\:\s+(.*)/ && do {
                    my $sub_ali = $1;
                    my $sub_sub = $2;
                    my %substituents = split(/[=\s]/, $sub_sub);
                    $sub->{$sub_ali} = { map { $_ - 1 => $substituents{$_} }
                                            keys %substituents };
                }
            }
        };
    }

    #combine ligand and sub information;
    for my $ligand (keys %{ $lig }) {
        $ligs_subs->{$ligand} = {};
        my @explicit_sub = grep { $_ =~ /^\d+$/ } 
                               keys %{ $lig->{$ligand}->{substituents} };
        my @inexplicit_sub = grep { $_ !~ /^\d+$/ }
                                keys %{ $lig->{$ligand}->{substituents} };
        my %explicit_sub = map { $_ => $lig->{$ligand}->{substituents}->{$_} }
                            @explicit_sub;
        my %inexplicit_sub = map { $_ => $lig->{$ligand}->{substituents}->{$_} }
                                @inexplicit_sub;

        #examine the inexplicit sub
        open (my $fh, "<$AARON/Subs/subs") or die "Cannot open AARON/Subs/subs";
        my %subs_record;

        while (<$fh>) {
            chomp;
            if ($_ =~ /[0-9a-zA-Z]/) {
                $subs_record{$_} = 1;    
            }
        }

        close($fh);

        for my $key (keys %inexplicit_sub) {
            unless (exists $subs_record{$key}) {
                print "The inexplicit substituent on the ligand $key " .
                      "Cannot found in our database " .
                      "This substituent is skipped.\n";
                delete $inexplicit_sub{$key};
            }
        }

        $ligs_subs->{$ligand}->{ligand} = $lig->{$ligand}->{ligand};
        $ligs_subs->{$ligand}->{inexp_sub} = \%inexplicit_sub;
        $ligs_subs->{$ligand}->{exp_sub} = \%explicit_sub;
        $ligs_subs->{$ligand}->{substrate} = {};
        $ligs_subs->{$ligand}->{jobs} = {};
        $ligs_subs->{$ligand}->{thermo} = {};

        for my $sub_ali (keys %{ $sub }) {
            my $lig_sub_ali = "$ligand-$sub_ali";
            $ligs_subs->{$ligand}->{substrate}->{$lig_sub_ali} = $sub->{$sub_ali};
        }
    }

    #chomp each value
    for my $key (keys %arg_in) {
        chomp $arg_in{$key};
    }
    close $in_h;
    
    #examine arguments;
    if ($#{$arg_in{selectivity}} == 0 &&
        ($arg_in{selectivity}->[0] eq 'NONE')) {
        $arg_in{selectivity} = [];
    }
    
    unless ($arg_in{template}) {
        print "A template must be figure out explicitly in the <jobname>.in file " .
              "by template=xxxxx, \n" .
              "If this catalyst contains different steps in a reaction, " .
              "you should figure out the step too. e.g. template=catalyst/TS1. \n" .
              "Exit without calculation\n";
        exit(1);
    }

    my $TS_path = (-d "$HOME/$TS_lib/$arg_in{reaction_type}/$arg_in{template}") ?
                    "$HOME/$TS_lib/$arg_in{reaction_type}/" : 
                    "$AARON/$TS_lib/$arg_in{reaction_type}/";

    unless (-d $TS_path) {
        print "Can't find $arg_in{template} in either user defined TS library: ".
              "$HOME/$TS_lib/ or the built_in library: $AARON/$TS_lib/$arg_in{template}\n";
        exit(1);
    }

    for my $level ( $arg_in{low_level}, $arg_in{level}, $arg_in{high_level} ) {
        $level->check_gen($arg_in{gen});
    }

    if ($arg_in{high_method}){
        unless ($arg_in{high_basis}->initiated()) {
            $arg_in{high_basis} = $arg_in{basis};
        }
    }

    $arg_in{TS_path} = $TS_path;
}


sub get_job_template {
    if ( -e "$AARON/template.job") {
        my $job_invalid;
        my $template_pattern = TEMPLATE_JOB;
        $template_job->{job} = "$AARON/template.job";
        $template_job->{formula} = {};
        open JOB, "<$AARON/template.job";
        #get formulas
        while (<JOB>) {
            /&formula&/ && do { while (<JOB>) {
                                    /&formula&/ && last;
                                    /^(\S+)=(\S+)$/ && do {  my $formula = $2;
                                                             my @pattern = grep {$formula =~ 
                                                                    /\Q$_\E/} values %$template_pattern;

                                                             unless (@pattern) {
                                                                print "template.job in $AARON is invalid. " .
                                                                      "Formula expression is wrong. " .
                                                                      "Please see manual.\n";
                                                                $job_invalid = 1;
                                                                last;
                                                            }
                                                            $template_job->{formula}->{$1} = $2 };
                                }
                                last if $job_invalid;
                              }
        }

    }
}


sub grab_cata_coords {
    #change catalyst from catalyst/TS into catalyst_TS
    my @temp = split('/', $arg_in{catalyst});
    my $catalyst = shift(@temp);
    $catalyst =~ s/\//_/g;
    my $cata_xyz = $catalyst.".xyz";
    unless (-e $cata_xyz) {
        print "You want to map new catalyst, but you didn't provide geometry of the new catalyst.\n";
        print "Exist without calculation.\n";
        exit(1);
    }
    my @cat_coords = grab_coords($cata_xyz);

    return @cat_coords;
}


sub write_status {
    my %status;

    my @ligands = keys %{ $ligs_subs };

    for my $ligand (@ligands) {
        my @geos = keys %{ $ligs_subs->{$ligand}->{jobs} };

        @status{@geos} = map {$ligs_subs->{$ligand}->{jobs}->{$_}->copy()} @geos;
    }

    my $string = Data::Dumper->Dump([\%status], ['STATUS']);

    open STATUS, ">.status" or die "Cannot open .status file: $!\n";

    print STATUS $string;

    close STATUS;
}


sub read_status {
    my $STATUS;

    my $string;
    if (-e ".status") {
        open STA, "<.status" or die "Cannot read status file: $!\n";
        $string = do {local $/; <STA>};

        eval $string;

        for my $key (keys %{ $STATUS }) {
            my ($head) = split(/\//, $key);
            ($head) = split(/\-/, $head);

            $ligs_subs->{$head}->{jobs}->{$key} = $STATUS->{$key};
        }
    }
}


#main function to initiate Aaron job
sub init_main {
    print "Preparing to run transition state searches...\n";
    sleep(2);
    &read_params();
    &get_job_template();
    &read_status();
    sleep(10);
}


package Theory_level;
use strict; use warnings;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

use AaronTools::Atoms qw(:BASIC);

sub new {
    my $class = shift;

    my %params = @_;

    my $self = { 
        basis => {},
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
        }else {
            if (@atoms) {
                @{$self->{basis}}{@atoms} = ($entry) x (scalar @atoms);
            }else {
                @{$self->{basis}}{keys %{ $masses }} = ($entry) x keys %{ $masses };
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

    my @gen_basis_atom = grep { $self->{basis}->{$_} =~ /gen\// } keys %{ $self->{atoms} };

    if (@gen_basis_atom && $gen) {
        for my $atom (@gen_basis_atom) {
            $self->{basis}->{$atom} =~ s/gen\///;
        }
        $self->{gen} = $gen;
    }elsif (@gen_basis_atom && (! $gen)) {
        print "You must provide path to the gen basis set if you want to use gen basis set\n";
        exit 0;
    }
}


sub method {
    my $self = shift;

    my $method;

    if ($self->{ecp}) {
        $method = $self->{method} . "/genecp";
    }elsif ($self->{gen}) {
        $method = $self->{method}. "/gen";
    }elsif (@{$self->unique_basis()} > 1) {
        $method = $self->{method}. "/gen";
    }elsif (@{$self->unique_basis()} == 1) {
        $method = $self->{method} . "/$self->unique_basis{}->[0]";
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

    if ($self->{gen}) {
        if (@{$self->unique_basis()} == 1) {
            my $basis = $self->unique_basis()->[0];
            $return = $self->{gen};
            $return =~ s/\/$//;
            $return .= "/$basis/N\n";
        }else {
            print "gen basis can only use one type of basis file. " .
                  "If you have multiple basis, please write them in the same gen basis file.";
            exit 0;
        }
    }elsif (@{$self->unique_basis()} > 1) {
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

            $return .= sprintf "%s " x @exit_atoms, @exit_atoms;
            $return .= "0\n";
            $return .= "$basis\n";
            $return .= '*' x 4 . "\n";
        }

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

            $return .= sprintf "%s " x @exit_atoms, @exit_atoms;
            $return .= "0\n";
            $return .= "$ecp\n";
        }
    }

    return $return;
}


1
