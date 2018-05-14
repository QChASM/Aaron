package AaronTools::JobControl;

use strict; use warnings;
use lib $ENV{'PERL_LIB'};
use lib $ENV{'AARON'};

use Constants qw(:SYSTEM :JOB_FILE);

use Exporter qw(import);
use Cwd qw(getcwd);

our @EXPORT = qw(findJob killJob submit_job count_time);

my $queue_type = $ENV{'QUEUE_TYPE'} ?
                 $ENV{'QUEUE_TYPE'} : ADA->{QUEUE};
my $g09root = $ENV{'G09_ROOT'} ?
              $ENV{'G09_ROOT'} : "/software/lms/g09_D01";       #absolute path to root directory for Gaussian09
                

#Returns jobIDs of all jobs (queued or running) in current directory, returns 0 if no jobs in queue match current directory
#This could be improved by searching for $Path more carefully!
#Works for PBS (default) or LSF
sub findJob {
    my $Path = $_[0];
    chomp($Path);

    #Strip leading directories off of $Path, to deal with different ways $HOME is treated
    $Path =~ s/^\S+$ENV{USER}//;
    
    my @jobIDs;

    if($queue_type =~ /LSF/i) {				#LSF queue
        my $bjobs=`bjobs -l 2> /dev/null`;
        #Combine into one line with no whitespace
        $bjobs =~ s/\s+//g;
        $bjobs =~ s/\r|\n//g;

        #First grab all jobs
        my @jobs = ($bjobs =~ m/(Job<\d+>.+?RUNLIMIT)/g);

        #parse each job looking for $Path
        foreach my $job (@jobs) {
            if ($job =~ /Job<(\d+)>\S+CWD<.+$Path>/) {
                push(@jobIDs,$1);
            }
        }
    }elsif ($queue_type =~ /PBS/i) {				#PBS (default)
        my $qstat = `qstat -fx`;
    
        #First grab all jobs
        my @jobs = ($qstat =~ m/<Job>(.+?)<\/Job>/g);
    
        #Grab jobIDs for all jobs matching $Path
        foreach my $job (@jobs) {
        	if ($job =~ m/<Job_Id>(\d+)\S+PBS_O_WORKDIR=\S+$Path</) {
        		push(@jobIDs, $1);
        	}
        }
    }elsif ($queue_type =~ /Slurm/i) {
        my @alljobs=`squeue -o %i_%Z`;
        foreach my $job (@alljobs) {
            if($job =~ /$Path/) {
                my @array = split(/_/, $job);
                push(@jobIDs,$array[0]);
            }
        }
    }
    
    if(@jobIDs) {
    	return @jobIDs;
    }
    return;
} #end sub findJob


sub killJob {
    my ($job) = @_;

    if ($queue_type =~ /LSF/i) {
        system("bkill $job");
    }elsif ($queue_type =~ /PBS/i) {
        system("qdel $job");
    }elsif ($queue_type =~ /Slurm/i) {
        system("scancel $job");
   }
   sleep(3);
} #end kill_job;


#Works for LSF or PBS (default)
sub submit_job {
    my %params = @_;

    my ($dir, $com_file, $walltime, 
        $numprocs, $template_job, $node) = ( $params{directory}, 
                                             $params{com_file},
                                             $params{walltime},
                                             $params{numprocs},
                                             $params{template_job},
                                             $params{node} );
    
    chomp(my $jobname=`basename $com_file .com`);
    my $jobfile = $dir ? "$dir/$jobname.job" : "$jobname.job";

    $dir //= getcwd();

    #if template_job were provide, use template_job to build job file
    if ($template_job->{job}) {

        my $job_file = $template_job->{job};

        my $job_found;

        my $template_pattern = TEMPLATE_JOB;

        open JOB_TEM, "<$job_file" or die "Cannot open $job_file:$!\n";
        my $job_content = do {local $/; <JOB_TEM>};
        close (JOB_TEM);

        $job_content =~ s/\Q$template_pattern->{JOB_NAME}/$jobname/g && do {$job_found = 1;};
        $job_content =~ s/\Q$template_pattern->{WALL_TIME}/$walltime/g;
        $job_content =~ s/\Q$template_pattern->{N_PROCS}/$numprocs/g;
        $job_content =~ s/\Q$template_pattern->{NODE_TYPE}/$node/g;
        #remove the formula part
        $job_content =~ s/&formula&\n(.*\n)*&formula&\n//g;

        for my $var (sort keys %{ $template_job->{formula} }) {
            my $var_value = eval($template_job->{formula}->{$var});
            $job_content =~ s/\Q$var/$var_value/g;
        }

        if ($job_found) {
            print "$jobfile\n";
            open JOB, ">$jobfile" or die "cannot open $jobfile\n";
            print JOB $job_content;
            close (JOB);
        }
    }

    
    unless (-e $jobfile) {

        my $memory=$numprocs*120000000;
        my $mb;
        if($queue_type eq 'LSF') {
            $memory = 0.8*$numprocs*2*10**9/8;		#memory in words per job
            $mb = 2700;				#memory in MB per core + padding
        }

        open JOB, ">$jobfile";

        if($queue_type =~ /LSF/i) {					#LSF queue
            print JOB "#BSUB -J $jobname\n" .
                      "#BSUB -o $jobname.job.%J\n" .
                      "#BSUB -L /bin/bash\n" .
                      "#BSUB -W $walltime:00\n" .
                      "#BSUB -M $mb\n" .
                      "#BSUB -R 'rusage[mem=$mb]'\n" .
                      "#BSUB -R 'span[ptile=$numprocs]'\n" .
                      "export g09root=$g09root\n" .
                      ". \$g09root/g09/bsd/g09.profile\n" .
                      "trap \"rm -r \$SCRATCH/\$LSB_JOBID\" 0 1 2 3 9 13 14 15\n" .
                      "mkdir \$SCRATCH/\$LSB_JOBID\n" .
                      "cd \$SCRATCH/\$LSB_JOBID\n" .
                      "echo -P- $numprocs > Default.Route\n" .
                      "echo -M- $memory >> Default.Route\n" .
                      "module purge\n" .
                      "env\n" .
                      "cp \$LS_SUBCWD/*.chk .\n" .
                      "g09  < \$LS_SUBCWD/$jobname.com  > \$LS_SUBCWD/$jobname.log\n" .
                      "cp *.chk \$LS_SUBCWD/\n" .
                      "exit\n";
        }elsif ($queue_type =~ /PBS/i) {					#PBS (default)
            print JOB "#!/bin/bash\n" .
                      "#PBS -l walltime=$walltime:00:00,mem=8gb,nodes=1:ppn=$numprocs\n\n\n" .
                      "export g09root=$g09root\n" .
                      ". \$g09root/g09/bsd/g09.profile\n\n" .
                      "trap \"\\rm -r \$TMPDIR/\$PBS_JOBID\" 0 1 2 3 9 13 14 15\n\n" .
                      "mkdir \$TMPDIR/\$PBS_JOBID\n" .
                      "cd \$TMPDIR/\$PBS_JOBID\n\n" .
                      "echo -P- $numprocs > Default.Route\n" .
                      "echo -M- $memory >> Default.Route\n\n" .
                      "module purge\n\n" .
                      "env\n\n" .
                      "cp \$PBS_O_WORKDIR/*.chk .\n" .
                      "g09  < \$PBS_O_WORKDIR/$jobname.com > \$PBS_O_WORKDIR/$jobname.log\n\n" .
                      "cp *.chk \$PBS_O_WORKDIR/\n\n" .
                      "exit";
        }elsif ($queue_type =~ /Slurm/i) {
            print JOB "#!/bin/bash\n" .
                       "#\n" .
                       "#SBATCH -J $jobname -e $jobname.job.e%j -o $jobname.job.o%j\n";
            if ($node) {
                print JOB "#SBATCH -p $node\n";
            }else {
                print JOB "#SBATCH -p medium\n";
            }
            print JOB "#SBATCH -t $walltime:00:00 -n $numprocs --mem=56G\n" .
                      "#SBATCH --ntasks-per-node=$numprocs\n" .
                      "\n" .
                      "\n" .
                      "cd \$TMPDIR\n" .
                      "\n" .
                      "df -h\n" .
                      "export g09root=/sw/group/lms/sw/g09_D01\n" .
                      ". $g09root/g09/bsd/g09.profile\n" .
                      "\n" .
                      "echo -P- 28 > Default.Route \n" .
                      "echo -M- 56GB >> Default.Route \n" .
                      "\n" .
                      "module purge \n" .
                      "\n" .
                      "env\n" .
                      "\n" .
                      "g09  <  \$SLURM_SUBMIT_DIR/$jobname.com  >  \$SLURM_SUBMIT_DIR/$jobname.log\n" .
                      "\n" .
                      "df -h\n" .
                      "\n" .
                      "ls -al\n" .
                      "\n" .
                      "\cp *.wf? *.47 \$LS_SUBCWD\n" .
                      "\n" .
                      "exit\n";
        }
        close(JOB);
    }

    my $failed = 1;
    #Alert user if qsub (or bsub) returns error
    #FIXME
    if (-e $jobfile) {
        my $current = getcwd();
        $failed = 0;

        chdir($dir);
        if($queue_type =~ /LSF/i) {
            if(system("bsub < $jobname.job >& /dev/null")) {
                print "Submission denied!\n";
                $failed = 1;
            }
        } elsif($queue_type =~ /Slurm/i) {
            if(system("qsub $jobname.job -N $jobname >& /dev/null")) { 
                print "Submission denied!\n";
                $failed = 1;
            }
        }
        chdir($current);
    }
    return $failed;
} #end sub submit


sub count_time {
    my ($sleep_time) = @_;

    my $time = localtime;
    
    my $sleep_hour = int($sleep_time/60);
    my $sleep_minute = $sleep_time%60;

    if ($time =~ /\s(\d+)\:(\d+)/) {
        my $hour = $1 + $sleep_hour;
        my $minute = $2 + $sleep_minute;
        
        if ($minute > 60) {
            $hour += 1;
            $minute += $minute%60;
        }

        $time =~ s/\s\d+\:\d+/ $hour:$minute/;
    }

    return $time;
}


sub call_g09 {
    my %params = @_;

    my ($com_file, $walltime, 
        $numprocs, $template_job,
        $node) = ($params{com_file}, $params{walltime}, $params{numprocs},
                  $params{template_job}, $params{node});

    chomp(my $jobname = `basename $com_file .com`);
    my $jobfile = "$jobname.job";

    my $shellfile = "$jobname.sh";

    unless (-e $shellfile) {

        open SHELL, ">$shellfile";
        
        my $template_pattern = TEMPLATE_JOB;

        my @job_command = @{$template_job->{command}};

        for my $command (@job_command) {
            $command =~ s/\Q$template_pattern->{JOB_NAME}/$jobname/g;
            $command =~ s/\Q$template_pattern->{WALL_TIME}/$walltime/g;
            $command =~ s/\Q$template_pattern->{N_PROCS}/$numprocs/g;
            $command =~ s/\Q$template_pattern->{NODE_TYPE}/$node/g;
            #remove the formula part
            $command =~ s/&formula&\n(.*\n)*&formula&\n//g;

            for my $var (sort keys %{ $template_job->{formula} }) {
                my $var_value = eval($template_job->{formula}->{$var});
                $command =~ s/\Q$var/$var_value/g;
            }

            next if ($command =~ /^exit/);
            print SHELL "$command\n";
        }
        close SHELL;

        chmod (0755, $shellfile);
    }

    my $walltime_sec = $walltime * 3600;
    eval {
        local $SIG{ALRM} = sub { die "TIMEOUT\n" };
        alarm $walltime_sec;
        eval {
            system("timeout $walltime_sec sh $shellfile");
        };
        alarm 0;
    };
    alarm 0;

    if ($@) {
        die unless $@ eq "TIMEOUT\n";
    }
}




  
    



1;
