package AaronInit;

use strict; use warnings;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

my $HOME = $ENV{'HOME'};
my $AARON = $ENV{'AARON'};

use Constants qw(:INFORMATION :THEORY :PHYSICAL :SYSTEM :JOB_FILE :OTHER_USEFUL);
use AaronTools::JobControl qw(get_job_template);
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

use Cwd qw(getcwd);
use Getopt::Long;
use Exporter qw(import);
use AaronTools::Atoms qw(:BASIC :LJ);
use AaronTools::G_Key;
use Data::Dumper;

#default values for some argument

our @EXPORT = qw($G_Key $W_Key %arg_parser $ligs_subs $parent $jobname
                 $template_job init_main write_status);

my $LIG_OLD = NAMES->{LIG_OLD};
my $LIG_NONE = NAMES->{LIG_NONE};
our $jobname;
our $parent = getcwd();
our $G_Key = new AaronTools::G_Key();
our $W_Key = new AaronTools::Workflow_Key();

my $input_file;

#content of template job file
our $template_job = get_job_template();

#ligand and substituent information
our $ligs_subs = {};

#read command line arguments
sub check_modules {
    print "Checking for required Perl modules...\n";

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
    
    print "Necessary Perl modules found. Starting AARON...\n"; 
}


        

sub read_args{
    GetOptions(
        'debug' => \$W_Key->{debug},
        'nosub' => \$W_Key->{nosub},
        'help|h' => \$arg_parser{help},
        'restart' => \$arg_parser{restart},
        'record' => \$W_Key->{record},
        'short' => \$W_Key->{short},
        'absthermo' => \$arg_parser{absthermo},
        'sleep=s' => \$arg_parser{sleeptime},
        'no_quota' => \$W_Key->{no_quota},
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

    ($jobname) = $input_file =~ /(\S+)\.in/;
}


#read arguments from input file
sub read_params {
   
    open my $in_h, "< $input_file" or die "Can't open $input_file:$!\n";

    ($jobname) = $input_file =~ /(\S+)\.in/; 

    my $lig = {};
    my $sub = {};
    my $custom;

    while(<$in_h>) {
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
    close $in_h;
    
    $G_Key->read_key_from_input($input_file);
    $W_Key->read_input($input_file);

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
        open (my $fh, "<$AARON/AaronTools/Subs/subs") or die "Cannot open $AARON/AaronTools/Subs/subs";
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
                print "The substituent on the ligand $key " .
                      "Cannot be found in our database " .
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
    
            $ligs_subs->{$head}->{jobs}->{$key}->{Gkey} = $G_Key;
            $ligs_subs->{$head}->{jobs}->{$key}->{Wkey} = $W_Key;
            $ligs_subs->{$head}->{jobs}->{$key}->{template_job} = $template_job;
            
            if ($ligs_subs->{$head}->{jobs}->{$key}->{conformers}) {
                for my $cf (sort keys %{ $ligs_subs->{$head}->{jobs}->{$key}->{conformers} }) {
                    $ligs_subs->{$head}->{jobs}->{$key}->{conformers}->{$cf}->{Gkey} = $G_Key;
                    $ligs_subs->{$head}->{jobs}->{$key}->{conformers}->{$cf}->{Wkey} = $W_Key;
                    $ligs_subs->{$head}->{jobs}->{$key}->{conformers}->{$cf}->{template_job} = $template_job;
                }
            }
        }
    }
}


#main function to initiate Aaron job
sub init_main {
    &read_args();
    &check_modules();
    print "Preparing to run transition state searches...\n";
    &read_params();
    &read_status();
    sleep(10);
}



1
