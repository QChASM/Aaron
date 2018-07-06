package Aaron::AaronOutput;

use strict;

use lib $ENV{'QCHASM'};

use Aaron::Constants qw(INFO);
use Aaron::AaronInit qw($G_Key %arg_parser $W_Key $parent $jobname $ligs_subs);

use Cwd;
use Exporter qw(import);

our @EXPORT = qw(&init_log print_message print_params terminate_AARON 
                 clean_up print_ee print_status close_log sleep_AARON);

my $QCHASM = $ENV{'QCHASM'};
my $ol;
my $out_file;
my $old_data;
my $queue = $ENV{'QUEUE_TYPE'};

sub init_log {
    my %params = @_;
    my ($job_name, $print_params) = ($params{job_name}, $params{print_params});
    $print_params //= 0;

    $job_name //= $jobname;

    $out_file = $parent . '/' . $job_name . ".log";
    if (-e $out_file) {
        open $ol, ">>$out_file" or die "Can't open $out_file\n";
        &restart_header($print_params);
    }else {
        open $ol, ">>$out_file" or die "Can't open $out_file\n";
        &header();
        &print_params() if $print_params;;
    }
}

#Prints Aaron header to $ol.  open STDOUT as $ol to write to screen
sub header {

    my $date = localtime;
    my $version = INFO->{VERSION};
    my @authors = @{ INFO->{AUTHORS} };
    my $year = INFO->{YEAR};

    &print_message("Aaron job started on $date\n\n");

    print $ol "                            Welcome to AARON!\n";
    print $ol "                                (v. $version)\n\n";
    print $ol "                               Written by\n";
    foreach my $author (@authors) {
      print $ol "                            $author\n";
    }
    print $ol "                          Texas A&M University\n";
    print $ol "                         September, 2013 - 2018\n\n";
    print $ol "                         University of Georgia\n";
    print $ol "                               2018 - \n\n";
    print $ol "                                          /\\          /\\ \n";
    print $ol "                                         ( \\\\        // )\n";
    print $ol "                                          \\ \\\\      // / \n";
    print $ol "                                           \\_\\\\||||//_/  \n";
    print $ol "                                            \\/ _  _ \\    \n";
    print $ol "                                           \\/|(O)(O)|    \n";
    print $ol "                                          \\/ |      |    \n";
    print $ol "                      ___________________\\/  \\      /    \n";
    print $ol "                     //  |          |  //     |____|     \n";
    print $ol "                    //   |A.A.R.O.N.| ||     /      \\    \n";
    print $ol "                   //|   |          | \\|     \\ 0  0 /    \n";
    print $ol "                  // \\   ||||||||||||  V    / \\____/     \n";
    print $ol "                 //   \\     /________(   _ /             \n";
    print $ol "                \"\"     \\   /    /    |  ||               \n";
    print $ol "                       /  /\\   /     |  ||               \n";
    print $ol "                      /  / /  /      \\  ||               \n";
    print $ol "                      | |  | |        | ||               \n";
    print $ol "                      |_|  |_|        |_||               \n";
    print $ol "                       \\_\\  \\_\\        \\_\\\\              \n\n\n";
    print $ol "                              Automated\n";
    print $ol "                              Alkylation\n";
    print $ol "                              Reaction\n";
    print $ol "                              Optimizer for\n";
    print $ol "                              New catalysts\n\n";
    print $ol "Citation:\n";
    print $ol "AARON, verson $version, Y. Guan, V. M. Ingman, B. J. Rooks, and S. E. Wheeler, Texas A&M University, $year.\n\n";
    print $ol "Y. Guan, V. M. Ingman, B. J. Rooks, and S. E. Wheeler, \"AARON: An Automated Reaction Optimizer for New Catlaysts\",\n J. Chem. Theory Comput. (submitted).\n\n";
    print $ol "The development of AARON is sponsored in part by the National Science Foundation,Grants CHE-1266022 and CHE-1665407.\n\n\n";
} #end sub header


sub restart_header {
    my $print_params = shift;
    my $date=localtime;
    if (-e $out_file) {
        print $ol "\n---------------------------------------------------------\nAaron job restarted on $date\n\n";
    } else {
        print $ol "Aaron job restarted on $date\n\n";
        &header();
        &print_params() if $print_params;
    }
}


sub print_message {
    print "$_[0]";
    print $ol "$_[0]";
}


#print all job parameters to $ol
sub print_params {
    my $version = INFO->{VERSION};
    my $AARON_HOME = "$QCHASM/Aaron";
    my $method = $G_Key->{level}->method();
    my $high_method = $G_Key->{high_level}->method();
    my $low_method = $G_Key->{low_level}->method();

    print $ol "----------------------------------------------------------------------------------\n";
    print $ol "Parameters\n";
    print $ol "----------------------------------------------------------------------------------\n";
    print $ol " AARON_HOME          = $AARON_HOME\n";
    print $ol "  version            = $version\n";
    print $ol "\n Reaction parameters:\n";
    print $ol "  reaction_type      = $W_Key->{reaction_type}\n" if $W_Key->{reaction_type};
    print $ol "  solvent            = $G_Key->{solvent}\n";
    print $ol "  temperature        = $G_Key->{temperature} K\n";
    print $ol "  MaxRTS             = $W_Key->{MaxRTS}\n" if $W_Key->{MaxRTS};
    print $ol "  MaxSTS             = $W_Key->{MaxSTS}\n" if $W_Key->{MaxSTS};
    print $ol "  TS_path            = $W_Key->{TS_path}\n" if $W_Key->{TS_path};
    print $ol "\n Methods:\n";
    print $ol "  method = $method\n";
    print $ol "  high level method  = $high_method\n" if $high_method;
    if(my $basis = $G_Key->{level}->footer_log()) {
        print $ol "  basis set file     = $basis\n";
    }
    print $ol "  solvent model      = $G_Key->{solvent_model}\n" if $G_Key->{solvent_model};
    print $ol "  low-level method   = $low_method\n";
    print $ol "\n Queue parameters:\n";
    print $ol "  wall               = $G_Key->{wall} hours\n";
    print $ol "  nprocs             = $G_Key->{n_procs}\n";
    print $ol "  shortwall          = $G_Key->{short_wall} hours\n" if $G_Key->{short_wall};
    print $ol "  shortprocs         = $G_Key->{short_procs}\n" if $G_Key->{short_procs};
    print $ol "  queue_name         = $queue\n" if $queue;

    if(@ARGV) {
        print $ol "\n command-line flags  = @ARGV\n";
    }
    print $ol "----------------------------------------------------------------------------------\n\n";
} #end sub print_params


sub print_status {

    print "\033[2J";                                              #clear the screen
    print "\033[0;0H";                                    #jump to 0,0
    my $date1=localtime;                          #current date
    
    my $running = 0;

    my $msg = "Status for all jobs...($date1)\n";


    for my $ligand (keys %{ $ligs_subs }) {

        my @start;
        my @running;
        my @done;
        my @finished;
        my @restart;
        my @pending;
        my @killed;
        my @skipped;
        my @repeated;
        my @sleeping;

        my $jobs = $ligs_subs->{$ligand}->{jobs};

        my $_get_job = sub {
            my $geo = shift;

            if (my ($cf) = $geo =~ /(Cf\d+)/) {
                $geo =~ s/\/Cf\d+//;
                return $jobs->{$geo}->{conformers}->{$cf};
            }else {
                return $jobs->{$geo}
            }
        };

        my $_print_status = sub {
            my $job = shift;

            my $geometry = $job->{name};

            STATUS: {
                if ($job->{status} eq 'start') { push(@start, $geometry);
                                                   last STATUS; }

                if ($job->{status} eq 'restart') { push(@restart, $geometry);
                                                     last STATUS; }

                if ($job->{status} eq 'pending') { push(@pending, $geometry);
                                                     last STATUS;}

                if ($job->{status} eq 'done') { push(@done, $geometry);
                                                  last STATUS;}

                if ($job->{status} eq 'finished') { push(@finished, $geometry);
                                                      last STATUS;}

                if ($job->{status} eq 'running') {push(@running, $geometry);
                                                    last STATUS;}

                if ($job->{status} eq 'killed') {push(@killed, $geometry);
                                                    last STATUS;}

                if ($job->{status} eq 'skipped') {push(@skipped, $geometry);
                                                        last STATUS;}

                if ($job->{status} eq 'repeated') {push(@repeated, $geometry);
                                                      last STATUS;}
                if ($job->{status} eq 'sleeping') {push(@sleeping, $geometry);
                                                       last STATUS;}
            }
        };

        $msg .= '-' x 80;
        $msg .= "\nStatus for $ligand jobs......\n";

        foreach my $geometry (sort keys %{ $jobs }) {
            if ($jobs->{$geometry}->{conformers}) {
                for my $cf (keys %{ $jobs->{$geometry}->{conformers} }) {
                    &$_print_status($jobs->{$geometry}->{conformers}->{$cf});
                }
            }else {
                &$_print_status($jobs->{$geometry});
            }
       }


        @start && do {$msg .= "\nThe following jobs are going to start:\n";};
        for my $geometry(@start) {
            $msg .= "$geometry is starting the AARON workflow using the geometry from the TS library\n";
        }

        @done && do {$msg .= "\nThe following jobs are done:\n";};
        for my $geometry(@done) {
            my $job = &$_get_job($geometry);
            my $step_done = $job->{step} - 1;
            $msg .= "$geometry step $step_done is done\n";
        }

        @finished && do {$msg .= "\nThe following AARON are finished: \n";};
        for my $geometry(@finished) {
            $msg .= "$geometry finished normally\n";
        }

        @running && do {$msg .= "\nThe following jobs are running:\n";};
        for my $geometry(@running) {
            my $job = &$_get_job($geometry);
            my $gradient = $job->{gout}->gradient();
            $msg .= "$geometry step $job->{step} attempt $job->{attempt} " .
                    "cycle $job->{cycle}. $gradient\n";
        }

        @pending && do {$msg .= "\nThe following jobs are pending:\n";};
        for my $geometry(@pending) {
            my $job = &$_get_job($geometry);
            $msg .= "$geometry step $job->{step} attempt $job->{attempt}: ";
            $msg .= ($job->{msg} or "No msg recorded") . "\n";
        }

        @restart && do {$msg .= "\nThe following jobs have been restarted for some reasons:\n";};
        for my $geometry(@restart) {
            my $job = &$_get_job($geometry);
            $msg .= "$geometry step $job->{step} " .
                    "restarted: ";
            if ($job->{msg}) {
                $msg .= $job->{msg};
            }else {
                $msg .= "No message recorded. ";
            }
            $msg .= "Now at attempt $job->{attempt}, cycle $job->{cycle}.\n";
        }

        @skipped && do {$msg .= "\nThe following jobs are skipped due to some error during the calculation: \n";};
        for my $geometry(@skipped) {
            my $job = &$_get_job($geometry);
            $msg .= "$geometry step $job->{step} " .
                    "was skipped by reason: \n";
            if ($job->{msg}) {
                $msg .= $job->{msg};
            }else {
                $msg .= "No skipped message recorded\n";
            }
        }

        @killed && do {$msg .= "\nThe following jobs are stopped:\n";};
        for my $geometry(@killed) {
            my $job = &$_get_job($geometry);
            $msg .= "$geometry\n step $job->{step} attemp $job->{attempt}: ";
            $msg .= ($job->{msg} or "No msg recorded") . "\n";
        }

        @repeated && do {$msg .= "\nThe following jobs are repeated conformers:\n";};
        for my $geometry(@repeated) {
            my $job = &$_get_job($geometry);
            $msg .= "$geometry $job->{msg}\n";
        }

        @sleeping && do {$msg .= "\nThe following jobs have not been started and are awaiting other jobs:\n";};
        for my $geometry(@sleeping) {
            $msg .= "$geometry\n";
        }

        print_message('=' x 80);
        print_message("\n$msg\n\n");

    #write status into .sta file

        if (@done || @running || @pending || @restart || @start) {
            $running = 1;
        }
    }
    return $running;
}


sub print_ee {
    my ($thermo, $absolute_only, $absolute) = @_;

    my $data = '';

    my $print_thermo = sub {
        my %params = @_;
        my ($geo, $thermo, $cf) = ($params{name}, $params{thermo}, $params{cf});


        if ($cf) {
            $data .= sprintf "%6s", $geo;
        }else {
            $data .= sprintf "%-6s", $geo;
        }

        foreach my $thermo (@{$thermo}) {
            if ($absolute_only || $arg_parser{multistep} || $absolute) {
                $data .= sprintf "%13.6f", $thermo;
            }else {
                $data .= sprintf "%10.1f", $thermo;
            }
        }
        $data .= "\n";
    };

    my @data_keys;
    if ($absolute_only) {
        @data_keys = sort grep {$thermo->{$_}->{found}} keys %{ $thermo };
    }else {
        @data_keys = sort grep {@{$thermo->{$_}->{sum}}} keys %{ $thermo };
    }

    unless (@data_keys) {
        $data .= "No data yet...\n";
        return $data;
    }else {
        $data .= "\n" . '~' x 90 . "\n";
    }



    my $ee = []; 
    my $er = {};
    if (keys %{ $thermo } == 2) {
        my ($sum1, $sum2) = map { $thermo->{$_}->{sum} } sort keys %{ $thermo };

        for my $i (0..7) {
            if ($sum1->[$i] && $sum2->[$i]) {
                $ee->[$i] = ($sum1->[$i] - $sum2->[$i])/($sum1->[$i] + $sum2->[$i]) * 100
            }
        }
    }elsif (keys %{ $thermo } > 2) {
        my $sum_total = [];

        my @sums = map { $thermo->{$_}->{sum} } sort keys %{ $thermo };

        for my $i (0..7) {
            for my $sum (@sums) {
                $sum_total->[$i] += $sum->[$i] if $sum->[$i];
            }
        }

        for my $key ( sort keys %{ $thermo } ) {
            $er->{$key} = [ map { $thermo->{$key}->{sum}->[$_]/$sum_total->[$_] }
                                (0..$#{ $thermo->{$key}->{sum} }) ];
        }
    }


    if ($G_Key->{high_level}->method()) {
        if ($absolute_only || $arg_parser{multistep} || $absolute) {
            $data .= sprintf "%19s%13s%13s%13s%13s%13s%13s%13s\n", 'E', 'H', 'G', 'G_Grimme',
                                                    'E\'', 'H\'', 'G\'', 'G_Grimme\'';
        }else {
            $data .= sprintf "%16s%10s%10s%10s%10s%10s%10s%10s\n", 'E', 'H', 'G', 'G_Grimme',
                                                    'E\'', 'H\'', 'G\'', 'G_Grimme\'';
        }
    }else {
        if ($absolute_only || $arg_parser{multistep} || $absolute) {
            $data .= sprintf "%19s%13s%13s%13s\n", 'E', 'H', 'G', 'G_Grimme',
        }else {
            $data .= sprintf "%16s%10s%10s%10s\n", 'E', 'H', 'G', 'G_Grimme',
        }
    }

    if (@$ee) {
        $data .= sprintf "%-6s", 'ee';
        for my $e (@$ee) {
            $data .= sprintf "%9.1f%%", $e;
        }
        $data .= "\n";
    }

    for my $key (@data_keys ) {
        $data .= sprintf "%-6s", $key;
        $data .= "\n";
        $data .= '-' x 86 . "\n";
        if (%{ $er }) {
            for my $e (@{ $er->{$key} }) {
                $data .= sprintf "%9.1f%%", $e;
            }
            $data .= "\n";
        }

        for my $geo (sort keys %{ $thermo->{$key}->{geos} }) {
            my $thermo_geo = $thermo->{$key}->{geos}->{$geo};

            $geo = (split(/\//, $geo))[-1];

            if ($absolute_only && $thermo_geo->{conformers}) {
                $data .= sprintf "%-6s\n", $geo;
            }else {
                &$print_thermo( name => $geo, 
                              thermo => $thermo_geo->{thermo}) if @{ $thermo_geo->{thermo} };
            }

            if ($thermo_geo->{conformers}) {
                for my $cf (sort keys %{ $thermo_geo->{conformers} }) {
                    &$print_thermo ( name => $cf,
                                   thermo => $thermo_geo->{conformers}->{$cf},
                                       cf => 1) if @{ $thermo_geo->{conformers}->{$cf} };
                }
            }
            $data .= '-' x 86 . "\n" if (!$absolute_only && @{ $thermo_geo->{thermo} });
        }
        $data .= "\n";
    }
    return $data;
}

sub clean_up {
    my ($workingfile) = @_;
    my $startdir = cwd;
    chdir ($workingfile) or die "Unable to enter dir $workingfile:$!\n";
    opendir(DIR, ".") or die "Unable to open dir $workingfile:$!\n";
    my @names = readdir(DIR) or die "Unable to read $workingfile:$!\n";
    closedir(DIR);
    for (@names) {
        /^\.+$/ && do {next;};

        -d && do {
            &clean_up($_);
            next;
        };

        /\.xyz/ && do {
            for (@names) {
                /\.job/ && do { system("rm -f $_"); next;};
            }
            next;
        }
    }
    chdir ($startdir);
} #End clean_up


sub close_log {
    my $date = localtime;
    print $ol "AARON stopped $date\n";
    close ($ol);
}


sub terminate_AARON {
    print_message("Aaron finished, terminated normally\n");
    clean_up($parent)
    &close_log();
}


sub sleep_AARON {
    print  "Aaron will check for next cycle after $arg_parser{sleeptime} seconds\n" .
           " Sleeping...\n";
    close ($ol);

    sleep($arg_parser{sleeptime});

    open $ol, ">>$out_file" or die "Can't open $out_file\n";
}

1;
