#!/usr/bin/perl -w

use strict;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

my $AARON = $ENV{'AARON'};

use G09Job;
use AaronInit qw($template_job);
use AaronTools::Geometry;
use AaronTools::G_Key;
use AaronOutput qw(init_log print_message close_log);

use Constants qw(:OTHER_USEFUL :SYSTEM :JOB_FILE :PHYSICAL);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $Gkey = new AaronTools::G_Key;
my $Wkey = new AaronTools::Workflow_Key();

my @basis;
my $method;
my @ecp;
my @chargemult;

my $restart;
my $sleep = SLEEP_TIME;
my $custom;
my $wall;
my $procs;

#read arguments
GetOptions(
    'restart' => \$restart,
    'record' => \$Wkey->{record},
    'method|m=s' => \$method,
    'wall|w=i' => \$wall,
    'n_procs|p=i' => \$procs,
    'node|n=i' => \$Gkey->{node},
    'basis|b=s' => \@basis,
    'ecp|e=s' => \@ecp,
    'sleep=s' => \$sleep,
    'denfit' => \$Gkey->{denfit},
    'chargemult|c=i{2}' => \@chargemult,
    'solvent|s=s' => \$Gkey->{solvent},
    'pcm=s' => \$Gkey->{pcm},
    'custom=s' => \$Gkey->{custom},
    'short' => \$Wkey->{short},
    'debug' => \$Wkey->{debug},
    'noquota' => \$Wkey->{no_quota},
    'emp_disp=s' => \$Gkey->{emp_disp},
);

my ($input_xyz) = grep { $_ =~ /\.xyz$/ } @ARGV;
my ($input_name) = $input_xyz =~ /(\S+)\.xyz/;

if (-e "$input_name.1.com") {
    $Gkey->read_key_from_com("$input_name.1.com");
    $Gkey->read_key_from_input();
}else {
    $Gkey->{charge} = shift @chargemult if $chargemult[0];
    $Gkey->{mult} = shift @chargemult if $chargemult[0];
    $Gkey->{level}->read_method($method);
    @basis && do{$Gkey->{level}->read_basis($_) for (@basis)};
    @ecp && do{$Gkey->{level}->read_ecp($_) for (@ecp)};
    #read uninitiated keywords from custom file
    $Gkey->read_key_from_input();
    #read gen
    $Gkey->{level}->check_gen($Gkey->{gen});
}    


my $geometry = new AaronTools::Geometry( name => $input_name );

#main 
init_log(job_name=>$input_name, print_params=>0);

my $wrong_con = $geometry->check_connectivity();

if (@$wrong_con) {
    my $msg = '';
    for my $atom (@$wrong_con) {
        $atom ++;
        $msg .= "Atoms too close around atom $atom, ";
    }
    $msg .= "Aaron has stopped this time, check your structure and fix. If this is right, ".
            "restart AARON with flag -nocrowd\n";
    print_message($msg);
    close_log();
    die "$msg";
}

my $G09job = new G09Job_TS_Single(
    name => $input_name,
    catalysis => $geometry,
    Gkey => $Gkey,
    Wkey => $Wkey,
    template_job => $template_job );

#unless (-e "$input_name.1.com") {
#    $G09job->build_com( directory => '.');
#    $G09job->set_status('2submit');
#}

while ($G09job->job_running()) {
    $G09job->check_status_run();
    $G09job->print_status();
    sleep($sleep);
}



