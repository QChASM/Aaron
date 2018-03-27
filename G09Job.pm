#Contributors: Yanfei Guan and Steven Wheeler

use lib $ENV{'PERL_LIB'};
use lib $ENV{'AARON'};

use Constants qw(:OTHER_USEFUL :PHYSICAL);

my $AARON = $ENV{'AARON'};
my $HART_TO_KCAL = HART_TO_KCAL;
my $MAXCYCLE = MAX_CYCLE;
my $launch_failed = 0;
my $CUTOFF = CUTOFF;


package G09Job;
use strict; use warnings;

use Cwd qw(cwd);
use File::Path qw(make_path);
use AaronInit qw(%arg_in %arg_parser $parent $system $template_job);
use AaronOutput qw(print_message close_logfile);
use AaronTools::G09Out;
use AaronTools::JobControl;
use Constants qw(:OTHER_USEFUL :COMPARE);
use Data::Dumper;


sub new {
    my $class = shift;
    my %params = @_;

    my ($name, $step, $cycle, $thermo,
        $attempt, $status, $catalysis) = ($params{name}, $params{step}, 
                                          $params{cycle}, $params{thermo},
                                          $params{attempt}, $params{status}, 
                                          $params{catalysis});
    $name //= '';
    $step //= 1;
    $cycle //= 1;
    $attempt //= 1;
    $status //= 'start';
    $catalysis //= {};
    $thermo //= [];

    my $self = {
        name => $name,
        step => $step,
        cycle => $cycle,
        attempt => $attempt,
        status => $status,
        catalysis => $catalysis,
        thermo => $thermo,
        msg => '',
        error => '',
    };

    bless $self, $class;

    return $self;
}


sub maxstep {
    my $self = shift;

    return $self->{maxstep};
}


sub copy {
    my $self = shift;

    my $new = $self->_copy();

    if ($self->{conformers}) {
        my @cfs = keys %{ $self->{conformers} };

        map { $new->{conformers}->{$_} = $self->{conformers}->{$_}->_copy() } @cfs;
    }

    return $new;
}


sub set_catalysis {
    my $self = shift;
    my ($catalysis) = @_;

    $self->{catalysis} = $catalysis;
} 


sub set_status {
    my $self = shift;
    my ($status) = @_;

    $self->{status} = $status;
}


sub set_msg {
    my $self = shift;
    my ($msg) = @_;

    $self->{msg} = $msg;
}


sub check_status_run {
    my $self = shift;


    if ($self->{conformers}) {
        my $lib = $self->{lib};

        #get status of each conformer
        for my $cf (sort keys %{ $self->{conformers} }) {
            if ($self->{conformers}->{$cf}->_skip_geometry()) {
                next;
            }
            $self->{conformers}->{$cf}->check_step();
        }

        #check for repeated conformers
        $self->check_conformers();
        $self->generate_conformers();
        #run calculations for each conformers
        for my $cf (sort keys %{ $self->{conformers} }) {
            if ($self->{conformers}->{$cf}->_skip_geometry(1)) {
                next;
            }
            $self->{conformers}->{$cf}->run_stepX();
        }
    }else {
        unless ($self->_skip_geometry(1)) {
            $self->check_step();
            $self->run_stepX();
        }
    }
}


sub _skip_geometry {
    my $self = shift;

    my ($check_finished) = @_; 

    my @status = ('repeated', 'killed');

    push (@status, 'finished') if $check_finished;

    if (grep { $_ eq $self->{status} } @status) {
        return 1;
    }else {
        return 0;
    }
}


sub _check_step {
    my $self = shift;

    my $maxstep = $self->maxstep();

    if ($self->{status} eq 'finished' && $self->{gout}) {return};

    my $geometry = $self->{name};

    my $file_name = "$geometry/" . $self->file_name();

    my $path = cwd;

    $path .= "/$geometry";

    my $jobrunning = AaronTools::JobControl::findJob($path);
    my $step = $maxstep;
    
    my $check_reaction;
    my $output;
    while($step >= 0) {

        if (-e "$file_name.$step.log") {
            $output = new AaronTools::G09Out( file => "$file_name.$step.log" );
            $self->{step} = $step;
            $self->{gout} = $output;
            my $geometry = $output->geometry();

            if ($output->finished_normal()) {
                $self->{status} = 'done';
                #update the catalysis geometry
                $self->{catalysis}->conformer_geometry($geometry);
                $check_reaction = 1;
            }elsif ($jobrunning) {
                $self->{status} = 'running';
                $self->{catalysis}->conformer_geometry($geometry);
                $check_reaction = 1;
            }elsif ($output->error()) {
                $self->{status} = 'failed';
                $self->{error} = $output->error();
                $self->{err_msg} = $output->{error_msg};
            }
            last;
        }elsif (-e "$file_name.$step.com") {
            $self->{step} = $step;
            if ($jobrunning) {
                $self->{status} = 'pending';
            }else {
                $self->{status} = '2submit';
            }
            last;
        }
        $step--;
    }
}


sub check_conformers {
    my $self = shift;
   
    my @cf = grep { $self->{conformers}->{$_}->{status} ne 'sleeping' &&
                    ($self->{conformers}->{$_}->{status} ne 'repeated') && 
                    ($self->{conformers}->{$_}->{status} ne 'pending') }
                  keys %{ $self->{conformers} };

    #sort @cf according to number
    map { $_ =~ s/Cf// } @cf; my @num = sort{ $a <=> $b } @cf;

    @cf = map { 'Cf' . $_ } @num;

    for my $i (0..$#cf-1 ) {
        if ($self->{conformers}->{$cf[$i]}->_cf_alive) {

            my $job_i = $self->{conformers}->{$cf[$i]};
            my $energy_i = $job_i->{gout}->{energy};
            my $geometry_i = $job_i->{gout}->{geometry};
            $job_i->{catalysis}->conformer_geometry($geometry_i);

            for my $j ($i+1..$#cf) {
                if ($self->{conformers}->{$cf[$j]}->_cf_alive) {

                    my $job_j = $self->{conformers}->{$cf[$j]};
                    my $energy_j = $job_j->{gout}->{energy};
                    my $geometry_j = $job_j->{gout}->{geometry};

                    if (abs($energy_i - $energy_j)*$HART_TO_KCAL < $CUTOFF->{E_CUTOFF}) {
                        my $rmsd = $job_i->{catalysis}->conformer_rmsd(ref_cata => $geometry_j);
                        if ($rmsd >= $CUTOFF->{RMSD_CUTOFF}) {next;}

                        #determine which one to kill
                        my $slow_geo = $j;

                        if ($job_i->{step} < $job_j->{step}) {
                            $slow_geo = $i;
                        }elsif ($job_i->{step} == $job_j->{step}) {
                            my $converged_i = grep { $job_i->{gout}->{gradient}->{$_}->{converged} eq 
                                                     'YES' } keys %{ $job_i->{gout}->{gradient} };
                            my $converged_j = grep { $job_j->{gout}->{gradient}->{$_}->{converged} eq 
                                                     'YES' } keys %{ $job_j->{gout}->{gradient} };
                            if ($converged_i < $converged_j) {
                                $slow_geo = $i;
                            }
                        }

                        if ($slow_geo == $i) {
                            $self->{conformers}->{$cf[$i]}->_repeated_cf($cf[$j]);
                            last;
                        }else {
                            $self->{conformers}->{$cf[$j]}->_repeated_cf($cf[$i]);
                        }

                    }
                }
            }
        }
    }
}


sub _cf_alive {

    my $self = shift;

    if ($self->{gout} &&
        $self->{gout}->{energy} &&
        $self->{gout}->{geometry} &&
        ($self->{status} ne 'repeated')) {
        return 1;
    }else {return 0;}
}


sub _repeated_cf {

    my $self = shift;

    my ($repeating_cf) = @_;

    $self->{status} = 'repeated';
    $self->{msg} = "repeat with $repeating_cf";

    $self->kill_running();

    delete $self->{gout};

    #system("rm -fr $self->{name}");
}


sub generate_conformers {
    my $self = shift;

    for my $cf (sort keys %{ $self->{conformers} }) {

        if ($self->{conformers}->{$cf}->_skip_geometry(1)) {
            next;
        }

        my $cf_obj = $self->{conformers}->{$cf};


        if ($cf_obj->{status} eq 'done' &&
            ($cf_obj->{step} >= 3)) {
            my $conf_geo = $cf_obj->{gout}->geometry()->copy();

            $cf_obj->{catalysis}->conformer_geometry($conf_geo);

            my $directory = $cf_obj->{name};

            my $file_name = "$directory/" . $cf_obj->file_name();

            my ($cf, $cf_num) = $directory =~ /(Cf(\d+))/;


            my $num_conf = $cf_obj->{catalysis}->number_of_conformers();
            my $conf_num = ($cf_num - 1) % $num_conf + 1;

            my $body_num = $cf_num - $conf_num;;
    
            my $array = $cf_obj->{catalysis}->conf_array( number => $conf_num );

            my @label_changed = grep { $array->[$_] > 1 } (0..$#{ $array });

            my $start_cf = $label_changed[-1] || 0;

            my $current_array = $array;

            for my $i ($start_cf..$#{ $array }) {
                if ($array->[$i] == 1) {
                    my $array_new = [@$array];
                    my $max_conf_number = $cf_obj->{catalysis}->max_conf_number($i);
                    for my $new_sub_conf_num (2..$max_conf_number) {
                        $array_new->[$i] = $new_sub_conf_num;

                        my $new_cf_num = $cf_obj->{catalysis}->conf_num($array_new);

                        $new_cf_num += $body_num;

                        my $new_cf = "Cf$new_cf_num";

                        if ($self->{conformers}->{$new_cf}->{status} eq 'sleeping') {
                            $cf_obj->{catalysis}->make_conformer( new_array => $array_new,
                                                                current_array => $current_array );
                            $cf_obj->{catalysis}->remove_clash();
                            $current_array = [@$array_new];

                            $cf_obj->{catalysis}->{conf_num} = $new_cf_num;
                            $self->{conformers}->{$new_cf}->{status} = '2submit';
                            $self->{conformers}->{$new_cf}->{step} = 2;
                            $self->{conformers}->{$new_cf}->build_com();
                            $self->{conformers}->{$new_cf}->{msg} = 'new conformer just started';
                        }
                    }
                }
            }
        }
    }
}


sub run_stepX {
    my $self = shift;

    my $directory = $self->{name};

    my $file_name = "$directory/" . $self->file_name();

    if ($self->{status} eq '2submit') {
        $self->submit();
        $self->{status} = 'pending';
    }elsif ($self->{status} eq 'done') {
        $self->{catalysis}->conformer_geometry($self->{gout}->{geometry});

        my $finished = $self->move_forward();

        unless ($finished) {
            $self->{step} ++;
            $self->{attempt} = 1;

            $self->build_com();

            $self->submit();
        }
    }elsif ($self->{status} eq 'failed') {
        if ($self->{gout}->{geometry}) {
            $self->{catalysis}->conformer_geometry($self->{gout}->{geometry});
        }else {
            my $com = "$file_name.$self->{step}.com";

            $self->{catalysis}->update_geometry($com);
        }
    
        $self->{attempt} ++;
        if ($self->{attempt} > 5) {
            $self->{status} = 'killed';
            $self->{msg} = 'killed because of too many attempts';
            return;
        }

        my $quit = $self->build_com();
        if ($quit) {
            chdir($parent);
            close_logfile();
            exit 0;
        }
        $self->submit();
    }
}


sub update_lib {
    my $self = shift;

    my $directory = $self->{name};

    my $filename = "$directory/" . $self->file_name();

    my ($head, $tail) = $self->_get_head_tail();

    my $lib_dir = $head . '_XYZ';

    $directory =~ s/$head/$lib_dir/;

    $directory =~ s/$tail//;

    if (! -d $directory) {
        make_path($directory);
    }

    my $repeated;
    if ($tail =~ /^Cf\d+$/) {
    #rearrange the name to start from 1
        opendir (DIR, $directory);
        my @files = readdir(DIR);
        @files = grep { $_ =~ /^Cf\d+\.xyz/ } @files;
        my $num = scalar @files + 1;
        $tail = "Cf" . $num;
        my $xyz = $self->{catalysis}->XYZ();
        @files = map { "$directory/" . $_ } @files;

        for my $file (@files) {
            open my $fh, "<$file" or die "Cannot open xyz file in library:$!\n";
            my $content = do { local $/; <$fh> };
            close $fh;
            if ($xyz eq $content) {
                $repeated = 1;
                last;
            }
        }
    }

    my $lib_file = $directory . "$tail" . ".xyz";

    if (! -e $lib_file && (! $repeated)) {
        $self->{catalysis}->printXYZ($lib_file);
    }
}


sub higher_level_thermo {
    my $self = shift;

    unless ($#{ $self->{thermo} } == 7) {
        for my $i (0..$#{ $self->{thermo} }) {
            $self->{thermo}->[$i+4] = $self->{gout}->energy() + 
                                      $self->{thermo}->[$i] - $self->{thermo}->[0];
        }
    }
}


sub _get_head_tail {
    my $self = shift;

    my $name = $self->{name};

    my @pattern = split(/\//, $name);

    return ($pattern[0], $pattern[-1]);
}


sub kill_running {
    my $self = shift;

    my $path = cwd;
    $path .= "/$self->{name}";

    for my $job(AaronTools::JobControl::findJob($path)) {
        AaronTools::JobControl::killJob($job);
    }
}


sub submit {
    if ($arg_parser{nosub}) {
        return;
    }

    my $self = shift;

    my $geometry = $self->{name};

    my $filename = "$geometry/" . $self->file_name();

    my $step = $self->{step};

    my ($wall, $nprocs);

    if ($step == 1) {
        $wall = $system->{SHORT_WALL};
        $nprocs = $system->{SHORT_PROCS};
    }else {
        $wall = $system->{WALL};
        $nprocs = $system->{N_PROCS};
    }

    if ($arg_parser{debug} || $arg_parser{short}) {
        $wall = $system->{SHORT_WALL};
    }

    unless($arg_parser{nosub}) {
        $launch_failed += AaronTools::JobControl::submit_job( 
                                 directory => $geometry,
                                  com_file => "$filename.$step.com",
                                  walltime => $wall,
                                  numprocs => $nprocs,
                              template_job => $template_job,
                                      node => $system->{NODE_TYPE} );
    }

    if ($launch_failed > MAX_LAUNCH_FAILED) {

        my $time = AaronTools::JobControl::count_time(60);

        my $msg = "AARON has failed to submit jobs to queue more than " .
                  MAX_LAUNCH_FAILED .
                  " times. AARON believe you may run out of hours or " .
                  "something wrong with the system. " .
                  "AARON will just sleep for one hour and continue.\n" .
                  "AARON will restart at $time\n" .
                  "sleeping...";

        print_message($msg);
        $launch_failed = 0;

        sleep(3600);
    }
}


sub file_name {
    my $self = shift;

    my $filename = $self->{name};

    $filename =~ s/\//\./g;

    return $filename;
}


sub build_com {
    my $self = shift;
    my %params = @_;

    my ($dir, $filename) = ($params{directory}, $params{name});

    $dir //= $self->{name};

    $filename //= $self->file_name();

    my $file_name = "$dir/$filename";

    my $low_method = $arg_in{low_level}->method();
    my $method = $arg_in{level}->method();
    my $high_method = $arg_in{high_level}->method();

    if ($arg_parser{debug}) {
        $method = $low_method;
        $high_method = $low_method;
    }

    if ($arg_in{denfit}) {
        $method .= ' denfit';
        $high_method .= ' denfit';
    }

    if ($arg_in{emp_dispersion}) {
        $method .= " EmpiricalDispersion=$arg_in{emp_dispersion}";
        $high_method .= " EmpiricalDispersion=$arg_in{emp_dispersion}";
    }

    my $step = $self->{step};
    my $error = $self->{error};
    my $catalysis = $self->{catalysis};

    my ($route, $footer, $print_flag) = $self->com_route_footer( 
        filename => $filename,
        low_method => $low_method,
        method => $method,
        high_method => $high_method,
    );

    ERROR: {    
        if ($error eq 'CONV') {my $scf_change = $route =~ /scf=xqc/ ?
                                                0 : ($route .= " scf=xqc");
            
                                my $message = " SCF convergence problems with $filename ";
                                if ($scf_change) {
                                    $message .= "...scf=xqc is in use\n";
                                }

                                $self->{msg} = $message;
                                $self->{status} = 'restart';
                                last ERROR;
                              }

        if ($error eq 'EIGEN') { $route =~ s/opt=\(/opt=\(noeigen,/;
                                 my $message = "Wrong number of negative eigenvalues for $filename. ";
                                 $message .= "...Adding noeigen.\n";
                                 $self->{msg} = $message;
                                 $self->{status} = 'restart';
                                 last ERROR;
                               }

        if ($error eq 'QUOTA') { if ( ! $arg_parser{no_quota}) {
                                    my $msg = "\nAARON thinks you hit the disk quota. "  .
                                           "Make more space or have quota increased. Then add " .
                                           "\"-noquota\" in the command line to restart\n";
                                    print_message($msg);
                                    return 1;
                                 }else{
                                     last ERROR;
                                 }
                               }

        if ($error eq "CHK" ||
            $error eq 'OPTSTEP') { open COM, "<$file_name.$step.com" or do { 
                                       my $msg = "$filename CHK problem, " .
                                                 "while Aaron cannot open $filename.$step.com. " .
                                                 "Aaron has skipped this, check that.\n";
                                       $self->{msg} = $msg;
                                       $self->{status} = 'skipped';
                                       close COM;
                                       last ERROR;
                                   };
                                   my $com_content = do {local $/; <COM>};
                                   close COM;
                                   my $msg = '';
                                   if (-e "$file_name.chk") {
                                       unlink "$file_name.chk";
                                       $msg = "Problem with check point file, using calcfc\n";
                                   }else {
                                       $msg = "No check point file was found in the directory, using calcfc\n";
                                   }
                                   $route =~ s/readfc/calcfc/;
                                   $self->{msg} = $msg;
                                   $self->{status} = 'restart';
                                   last ERROR;
                                 }
        if ($error eq "CLASH") { my $msg = "Atoms too crowded in $filename.$step.com. ";
                                 unless ($catalysis->remove_clash()) {
                                    $msg .= "Aaron failed to remove clash, please try manually\n";
                                    $self->{msg} = $msg;
                                    $self->{status} = 'skipped';
                                    
                                    return 0;
                                 }else {
                                     $msg .= "Aaron has removed the clash, the job is restarted. ";
                                     $self->{msg} = $msg;
                                     $self->{status} = 'restart';
                                 }
                               }
        if ($error eq "CHARGEMULT") { my $msg = "The combination of multipicity is not correct " .
                                                "AARON believe this is a fatal error so AARON quit " .
                                                "at this point, after you fix the problem, restart AARON\n";
                                      print_message($msg);
                                      return 1;
                                    }

        if ($error eq "REDUND") { my $msg = "Bend failed for some angles for $filename, " .
                                            "Aaron will restart this job.\n";

                                  $self->{msg} = $msg;
                                  $self->{status} = 'restart';
                                  if ($arg_parser{record}) {
                                      system("mv $file_name.$step.log $file_name.log.$step");
                                  }else{
                                      system("rm -fr $file_name.$step.*");
                                  }
                                }

        if ($error eq "UNKNOWN") { my $msg = "unknown reason, " .
                                             "AARON will retry the failed step. ";
                                   $msg .= "Please also check $filename.$step.log manually\n";

                                   $self->{msg} = $msg;
                                   $self->{status} = 'restart';
                                   if ($arg_parser{record}) {
                                       system("mv $file_name.$step.log $file_name.log.$step");
                                   }else{
                                       system("rm -fr $file_name.$step.*");
                                   }
                                 }

    }

    if ($self->{cycle} > 1) {
        my $fc_change = ($route =~ s/readfc/calcfc/);
        ($route, my $step_change) = &reduce_maxstep($route);
        if ($fc_change) {
            $self->{msg} .= "calculate fc instead of read fc from .chk file.\n";
            system("rm -fr $file_name.chk");
        }
        if ($step_change) {
            $self->{msg} .= "using smaller step.\n";
        }
    }

    if ($self->{attempt} > 2) {
        my $fc_change = ($route =~ s/readfc/calcfc/);
        ($route, my $step_change) = &reduce_maxstep($route);
        if ($step_change) {
            $self->{msg} .= "using smaller step.\n";
        }
        if ($self->{attempt} > 3) {
            unless ($route =~ /nonlinear/) {$route =~ s/opt=\(/opt=\(nolinear,/;}
            if ($fc_change) {
                $self->{msg} .= "calculate fc instead of read fc from .chk file.\n";
                system("rm -fr $file_name.chk");
            }
            if ($self->{attempt} > 4) {
                my $fc_change = (($route =~ s/readfc/calcfcall/) || ($route =~ s/calcfc/calcfcall/));
                if ($fc_change) {
                    $self->{msg} .= "calculate fc before each optimization step, ".
                                    "this can take long time.\n";
                    system("rm -fr $file_name.chk");
                }
            }
        }
    }

    my $com_file = "$file_name.$step.com";
    my $comment = "step $step (attempt $self->{attempt}) cycle $self->{cycle}";
    $self->{catalysis}->write_com( comment => $comment,
                                    route => $route,
                                    charge => $arg_in{charge},
                                    mult => $arg_in{mult},
                                    footer => $footer,
                                    print_flag => $print_flag,
                                    filename => $com_file );
    return 0;
} #End build new com file



#subroutine reduce maxstep from maxstep 
sub reduce_maxstep {
    my ($route_temp) = @_ ;
    my $step_change;
    if ($route_temp =~ /maxstep=(\d+)/) {
        my $new_max_step = $1 - 2;
        if ($new_max_step > 0) {
            $route_temp =~ s/maxstep=$1/maxstep=$new_max_step/;
            $step_change = 1;
        }
    }else {
        $route_temp =~ s/opt=\(/opt=\(maxstep=5,/;
        $step_change = 1;
    }

    return ($route_temp, $step_change);
}


package G09Job_TS;
use strict; use warnings;

use Cwd qw(cwd);
use AaronInit qw(%arg_in %arg_parser $parent $system $template_job);
use AaronOutput qw(print_message close_logfile);
use AaronTools::G09Out;
use AaronTools::JobControl;
use Constants qw(:OTHER_USEFUL);

our @ISA = qw(G09Job);

sub new {
    my $class = shift;
     
    my $self = new G09Job(@_);

    bless $self, $class;

    $self->{maxstep} = MAXSTEP->{TS};
    $self->{maxstep}-- unless $arg_in{high_method};

    return $self;
}


sub new_conformer {
    my $self = shift;
    my %params = @_;

    my $cf = $params{conf};

    $self->{conformers} = {} unless $self->{conformers};

    $self->{conformers}->{$cf} = new G09Job_TS( 
              name => $self->{name} . "/$cf",
              step => $params{step},
             cycle => $params{cycle},
           attempt => $params{attempt},
         catalysis => $self->{catalysis}
    );

    $self->{conformers}->{$cf}->set_status('sleeping');
}


#This is used to write status in the .status file in the AaronInit module
#NOTE this copy method will not copy the catlysis attribute
sub _copy {
    my $self = shift;

    my $new = new G09Job_TS( name => $self->{name},
                             step => $self->{step},
                            cycle => $self->{cycle},
                           attemp => $self->{attempt},
                           status => $self->{status},
                              msg => $self->{msg},
                           thermo => [@{ $self->{thermo} }],
                            error => $self->{error} );
    return $new;
}


sub check_step {
    my $self = shift;

    $self->_check_step();

    if ($self->{status} eq 'running' ||
        ($self->{status} eq 'done')) {
        $self->check_reaction() if ($self->{step} > 2);
    }
}


sub check_reaction {

    my $self = shift;

    my $geometry = $self->{name};

    my $filename = "$geometry/" . $self->file_name();

    my $catalysis = $self->{catalysis};

    my ($failed, $con) = $catalysis->examine_constraints();

    if ($failed) {
        $self->kill_running();

        $self->{cycle} ++;

        if ($self->{cycle} > $MAXCYCLE) {
            $self->{status} = 'killed';
            $self->{msg} = "killed because of too many cycles. ";
            return;
        }

        $self->{catalysis}->update_geometry("$filename.2.log");

        $self->remove_later_than2(); 

        $self->{step} = 2;
        $self->{attempt} = 1;

        my $distance = $failed * 0.1;

        $self->{catalysis}->change_distance( atom1 => $con->[0],
                                             atom2 => $con->[1],
                                       by_distance => $distance );
        $self->{catalysis}->_update_geometry();

        print_message("Changing the distance by $distance A\n");

        $self->build_com();

        $self->{status} = '2submit';
        $self->{msg} = "reverted to step 2, now waiting in the queue ";
    }
}


sub move_forward {
    my $self = shift;

    my $finished = 0;

    if ($self->{step} == 3) {
        $self->update_lib();
    }

    if ($self->{step} >= 4) {
        $self->get_thermo();
        if ($self->{step} == $self->maxstep()) {
            $self->higher_level_thermo() if $arg_in{high_method};
            $self->{status} = 'finished';
            $finished = 1;
        }
    }

    return $finished;
}


sub get_thermo {
    my $self = shift;

    my $directory = $self->{name};

    my $filename = "$directory/" . $self->file_name();

    my $out = $self->{gout};

    if (! @{ $self->{thermo} }) {
        if (! $out->{enthalpy}) {
            $out = new AaronTools::G09Out( file => "$filename.4.log" );
        }
        $self->{thermo} = [ $out->energy(), $out->enthalpy(),
                            $out->free_energy(), $out->Grimme_G() ];
    }
}


sub remove_later_than2 {

    my $self = shift;

    my $geometry = $self->{name};

    my $filename = "$geometry/" . $self->file_name();

    #remove evering thing more than step 2
    foreach my $later_step (2..$self->maxstep()) {
        if (-e "$filename.$later_step.com") {
            print_message("Removing $filename.$later_step.com...\n");
            if ($arg_parser{record}) {
                system("mv $filename.$later_step.com $filename.com.$later_step");
            }else {
                system("rm -fr $filename.$later_step.com");
            }
        }

        if (-e "$filename.$later_step.log") {
            print_message("Removing $filename.$later_step.log...\n");
            if ($arg_parser{record}) {
                system("mv $filename.$later_step.log $filename.$later_step.log.$later_step");
            }else {
                system("rm -fr $filename.$later_step.log");
            }
        }

        if (-e "$filename.$later_step.job") {
            print_message("Removing .job files for $later_step step...\n");
            system("rm -fr $filename.$later_step.job*");
        }
    }
}


sub com_route_footer {

    my $self = shift;

    my %params = @_;

    my ($low_method, $method, 
        $high_method, $filename) = ( $params{low_method},
                                     $params{method}, 
                                     $params{high_method}, 
                                     $params{filename} );

    my $dir //= $self->{name};

    $filename //= $self->file_name();

    my $file_name = "$dir/$filename";

    my $step = $self->{step};

    my $route;
    my $footer;
    my $catalysis = $self->{catalysis};

    my $print_flag;

    SWITCH: {
        if ($step == 1) { $route .= "#$low_method opt nosym";
                          #add constrats to substrate and old part of catalyst
                          $print_flag = 1;
                          $footer = $arg_in{low_level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 2) { $route .= "#$method opt=(modredundant,maxcyc=1000)";
                            
                          for my $constraint (@{ $catalysis->{constraints} }) {
                              my @bond = map { $_ + 1 } @{$constraint->[0]};
                              $footer .= "B $bond[0] $bond[1] F\n";
                          }
                          $footer .= "\n";
                          $footer .= $arg_in{level}->footer($catalysis);

                          last SWITCH; }

        if ($step == 3) { $route = "\%chk=$filename.chk\n";
                          if (-f "$file_name.chk") {
                            $route .= "#$method opt=(readfc,ts,maxcyc=1000)";
                          }else {
                            $route .= "#$method opt=(calcfc,ts,maxcyc=1000)";
                          }
                          $footer .= $arg_in{level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 4) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$method freq=(hpmodes,noraman,temperature=$arg_in{temperature})";
                          $footer .= $arg_in{level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 5) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$high_method";
                          $footer .= $arg_in{high_level}->footer($catalysis);
                          last SWITCH; }
    }

    if ($arg_in{solvent} ne "gas" && $step > 1) {
        $route .= " scrf=($arg_in{pcm},solvent=$arg_in{solvent})";
    }

    return ($route, $footer, $print_flag);
}


package G09Job_MIN;
use strict; use warnings;

use Cwd qw(cwd);
use AaronInit qw(%arg_in %arg_parser $parent $system $template_job);
use AaronOutput qw(print_message close_logfile);
use AaronTools::G09Out;
use AaronTools::JobControl;
use Constants qw(:OTHER_USEFUL);

our @ISA = qw(G09Job);

sub new {
    my $class = shift;

    my $self = new G09Job(@_);

    bless $self, $class;

    $self->{maxstep} = MAXSTEP->{INT};
    $self->{maxstep}-- unless $arg_in{high_method};

    return $self;
}


sub new_conformer {
    my $self = shift;
    my %params = @_;

    my $cf = $params{conf};

    $self->{conformers} = {} unless $self->{conformers};

    $self->{conformers}->{$cf} = new G09Job_MIN( 
              name => $self->{name} . "/$cf",
              step => $params{step},
             cycle => $params{cycle},
           attempt => $params{attempt},
         catalysis => $self->{catalysis}
    );

    $self->{conformers}->{$cf}->set_status('sleeping');
}


#This is used to write status in the .status file in the AaronInit module
#NOTE this copy method will not copy the catlysis attribute
sub _copy {
    my $self = shift;

    my $new = new G09Job_MIN( name => $self->{name},
                              step => $self->{step},
                             cycle => $self->{cycle},
                            attemp => $self->{attempt},
                            status => $self->{status},
                               msg => $self->{msg},
                            thermo => [@{ $self->{thermo} }],
                             error => $self->{error} );
    return $new;
}


sub check_step {
    my $self = shift;

    $self->_check_step();
}


sub move_forward {
    my $self = shift;

    my $finished = 0;

    if ($self->{step} == 2) {
        $self->update_lib();
    }

    if ($self->{step} >= 3) {
        $self->get_thermo();
        if ($self->{step} == $self->maxstep()) {
            $self->higher_level_thermo() if $arg_in{high_method};
            $self->{status} = 'finished';
            $finished = 1;
        }
    }

    return $finished;
}


sub get_thermo {
    my $self = shift;

    my $directory = $self->{name};

    my $filename = "$directory/" . $self->file_name();

    my $out = $self->{gout};

    if (! @{ $self->{thermo} }) {
        if (! $out->{enthalpy}) {
            $out = new AaronTools::G09Out( file => "$filename.3.log" );
        }
        $self->{thermo} = [ $out->energy(), $out->enthalpy(),
                            $out->free_energy(), $out->Grimme_G() ];
    }
}


sub com_route_footer {

    my $self = shift;

    my %params = @_;

    my ($low_method, $method, 
        $high_method, $filename) = ( $params{low_method},
                                     $params{method}, 
                                     $params{high_method}, 
                                     $params{filename} );

    my $dir //= $self->{name};

    $filename //= $self->file_name();

    my $file_name = "$dir/$filename";

    my $step = $self->{step};

    my $footer;
    my $route;

    my $catalysis = $self->{catalysis};

    my $print_flag;

    SWITCH: {
        if ($step == 1) { $route .= "#$low_method opt nosym";
                          #add constrats to substrate and old part of catalyst
                          $print_flag = 1;
                          $footer = $arg_in{low_level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 2) { $route = "\%chk=$filename.chk\n";
                          if (-f "$file_name.chk") {
                            $route .= "#$method opt=(readfc,maxcyc=1000)";
                          }else {
                            $route .= "#$method opt=(calcfc,maxcyc=1000)";
                          }
                          $footer .= $arg_in{level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 3) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$method freq=(hpmodes,noraman,temperature=$arg_in{temperature})";
                          $footer .= $arg_in{level}->footer($catalysis);
                          last SWITCH; }

        if ($step == 4) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$high_method";
                          $footer .= $arg_in{high_level}->footer($catalysis);
                          last SWITCH; }
    }

    if ($arg_in{solvent} ne "gas" && $step > 1) {
        $route .= " scrf=($arg_in{pcm},solvent=$arg_in{solvent})";
    }

    return ($route, $footer, $print_flag);
}










