#!/usr/bin/perl -w
use strict; use warnings;

my $failed_to_submit;
my $job_found;
my $failed_to_kill;

eval {
    use lib $ENV{'AARON'};
    use lib $ENV{'PERL_LIB'};
   
    use AaronTools::G_Key;
    use AaronTools::JobControl qw(get_job_template submit_job findJob killJob);

    my $Gkey = new AaronTools::G_Key;

    $Gkey->read_key_from_input();

    my $wall = $Gkey->{wall};
    my $n_procs = $Gkey->{n_procs};
    my $node = $Gkey->{node};

    my $template_job = get_job_template();

    $failed_to_submit = submit_job(
        com_file=>'test.com',
        walltime=>$wall,
        numprocs=>$n_procs,
        template_job=>$template_job,
        node=>$node);

    sleep(10);

    ($job_found) = findJob('.');

    killJob($job_found) if $job_found;
    sleep(5);

    $failed_to_kill = findJob('.');
    1
} or do {
    my $error = $@;

    die "Error found in code: $error\n";
};

if ($failed_to_submit) {
    die "Failed to submit test job to the queue.\n"
}

unless ($job_found) {
    die "Cannot find job submitted to the queue.\n";
}

if ($failed_to_kill) {
    die "Cannot kill the job on the queue.\n";
}

system("rm -fr test.job*");
system("rm -fr test.log");

print "Test passed!\n";



