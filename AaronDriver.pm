package Aaron::AaronDriver;

use strict;

use lib $ENV{'QCHASM'};

use AaronTools::Constants qw(PHYSICAL NAMES UNIT);
use Aaron::AaronInit qw($G_Key %arg_parser $parent
                 $ligs_subs $template_job $W_Key);

use Aaron::AaronOutput qw(print_message print_message_to_log print_ee print_to_thermo);
use AaronTools::Catalysis;
use Aaron::G09Job;

use File::Path qw(make_path);
use Cwd qw(getcwd cwd);
use Exporter qw(import);
use Math::Trig;
use File::Path qw(make_path);
use File::Copy qw(cp);
use Data::Dumper;

our @EXPORT = qw(make_directories get_status_run print_status analyze_result);


my $sub_ref;
my $constraints;
my @cat_coords;
my $thermo = {};
my $cata_read = {};

my $launch_failed = 0;

my $LIG_OLD = NAMES->{LIG_OLD};
my $LIG_NONE = NAMES->{LIG_NONE};
my $SUBSTRATE = NAMES->{SUBSTRATE};
my $LIGAND = NAMES->{LIGAND};

my $hart_to_kcal = UNIT->{HART_TO_KCAL};

my $queue_type = $ENV{'QUEUE_TYPE'};
my $skip_step1;

my $ts_found;
my $min_found;

#################
#ref to some useful function
################
my $imitate;
$imitate = sub {
    my ($working_dir, $f_file, $f_dir) = @_;
    my $start_dir = cwd;
    chdir($working_dir) or die "Unable to enter dir $working_dir: $!\n";
    opendir(DIR, '.') or die "Unable to open dir $working_dir: $!\n";
    my @names = sort readdir(DIR) or die "Unable to read $working_dir: $!\n";
    closedir(DIR);

    my %step_names = map{ $_ => 1 } @{ $W_Key->{step} };
    my @names_in_step = grep { exists $step_names{$_} } @names;
    @names_in_step && do { @names = @names_in_step };

    foreach my $name (@names) {
        next if ($name eq '.');
        next if ($name eq '..');

        if (-d $name) {
            &$imitate($name, $f_file, $f_dir);
            next;
        }
        &$f_file($name);
    }
    chdir($start_dir);
};  #end of $imitate


sub make_directories {
    for my $ligand(keys %{ $ligs_subs }) {
        &_make_directories($ligand);
    }
}


#initiate a a directory tree and status hashes
sub _make_directories {

    my ($lig_ali) = @_;

    if (!-d $lig_ali) {
        mkdir "$lig_ali";
    }
    chdir($lig_ali) or die "Unable to change into $lig_ali: $!\n";

    if ($ligs_subs->{$lig_ali}->{ligand} eq $LIG_NONE) {
        for my $sub (keys %{ $ligs_subs->{$lig_ali}->{substrate} }) {
            if (!(-d $sub)) {
                mkdir "$sub";
            }

            &dir_tree( target => $W_Key->{TS_path} . "$W_Key->{template}",
                       substrate => $sub,
                       ligand => $lig_ali,
                  no_new_subs => $W_Key->{input_conformers_only} );
        }

        unless (%{ $ligs_subs->{$lig_ali}->{substrate} }) {
            if ($W_Key->{input_conformers_only}) {
                my $msg = "Aaron will reoptimize template structures using the specified level of theory\n";
                print_message($msg);
                print_message_to_log($msg);
                $skip_step1 = 1;
            }else{
                my $msg = "Initiating conformer searching for $lig_ali (original catalyst).\n";
                print_message($msg);
                print_message_to_log($msg);
            }
            if (!(-d $lig_ali)) {
                mkdir "$lig_ali";
            }
            &dir_tree( target => $W_Key->{TS_path} . "$W_Key->{template}",
                       ligand => $lig_ali,
                  no_new_subs => $W_Key->{input_conformers_only} );
        }
    }else {
        if (!-d $lig_ali) {
            mkdir "$lig_ali";
        }
        &dir_tree( target => $W_Key->{TS_path} . "$W_Key->{template}",
                   ligand => $lig_ali,
              no_new_subs => $W_Key->{input_conformers_only} );

        for my $sub (keys %{ $ligs_subs->{$lig_ali}->{substrate} }) {
            if (!(-d $sub)) {
                mkdir "$sub";
            }

            my $top_dir_tar = "$parent/$lig_ali/$lig_ali" . "_XYZ";
            &dir_tree( target => $top_dir_tar,
                       substrate => $sub,
                       ligand => $lig_ali,
                   no_new_subs => 1);
        }
    }

    chdir($parent);
}


#mapping a directory tree and do any operation related to the specific file
sub dir_tree {
    my $new_dir = {};
    #top directory to mimic; top directory to make;
    #constraints for each gemotry;
    #status for each geometry;
    my %params = @_;
    my ($top_dir_tar, $substrate,
        $ligand, $no_new_subs) = ( $params{target},
                                   $params{substrate},
                                   $params{ligand},
                                   $params{no_new_subs} );
    if (!(-d $top_dir_tar)) {
        return;
    }

    my $substituents = {};
    my $substituents_ligand;
    my $inexp_sub;
    my $top_dir_make;
    my $new_ligand;
    my $status = $ligs_subs->{$ligand}->{jobs};

    my $submission_template = $template_job;
    #use the default $QCHASM/AaronTools/template.job unless a submission template was specified in the input file
    if( $W_Key->{submission_template} ) {
        $submission_template = $W_Key->{submission_template_o};
    }

    if ($substrate) {
        if (! -d $substrate) {
            mkdir $substrate;
        }
        $top_dir_make = "$substrate";
        $substituents->{substrate} =
            $ligs_subs->{$ligand}->{substrate}->{$substrate};
    }else {
        if (! -d $ligand) {
            mkdir $ligand;
        }
        $top_dir_make = "$ligand";
        if ($ligs_subs->{$ligand}->{ligand} ne $LIG_OLD &&
            ($ligs_subs->{$ligand}->{ligand} ne $LIG_NONE)) {
            $new_ligand = $ligs_subs->{$ligand}->{ligand};
            $substituents_ligand =
                $ligs_subs->{$ligand}->{exp_sub};
        }else {
            $substituents->{ligand} =
                $ligs_subs->{$ligand}->{exp_sub};
        }
        $inexp_sub = $ligs_subs->{$ligand}->{inexp_sub};
    }

    my $current_dir = cwd;

    chdir($top_dir_tar) or die "Unable to enter dir $top_dir_tar: $!\n";

    my $f_file = sub {
        my ($name) = @_;
        if ($name =~ /(\S+).xyz/) {
            my $extend = $1;
            my $tempdir = cwd;

            $ts_found = 1 if ($extend =~ /^ts/i);
            $min_found = 1 if ($extend =~ /^min/i);

            if ($tempdir =~ /$top_dir_tar(\S+)?/) {

                my $newdir = $1 ? $top_dir_make . $1 . "/$extend" : $top_dir_make . "/$extend";
                my $head = $newdir; $head =~ s/\/Cf\d+$//;

                unless ($cata_read->{$newdir}) {
                    print "  Checking $newdir...\n";
                    #make distance hashes for each geometry
                    my $catalysis = new AaronTools::Catalysis( name => $extend,
                                                       substituents => $substituents,
                                                        no_new_subs => $no_new_subs );

                    $catalysis->freeze_all_atoms();
                    if ($new_ligand) {
                        my $ligand = new AaronTools::Ligand( name => $new_ligand,
                                                     substituents => $substituents_ligand,
                                                      no_new_subs => $no_new_subs );
                        $catalysis->map_ligand($ligand);
                    }

                    $catalysis->substitute();
                    for my $sub (keys %{ $inexp_sub }) {
                        $catalysis->substitute( component => $LIGAND,
                                                target => $sub,
                                                sub => $inexp_sub->{$sub} );
                    }

                    $new_dir->{$newdir}->{catalysis} = $catalysis;
                    $new_dir->{$newdir}->{extend} = $extend;
                    $cata_read->{$newdir} = 1;
                }
            }
        }
    };

    &$imitate($top_dir_tar, $f_file);

    $W_Key->{multistep} = 1 if ($ts_found && $min_found);

    chdir($current_dir);

    #make paths and initialize every geometry to be on step0 1st attempt
    foreach my $newdir (keys %{ $new_dir }) {
        my $cat_temp = $new_dir->{$newdir}->{catalysis};

        my $cf_num = $cat_temp->number_of_conformers();

        my $extend = $new_dir->{$newdir}->{extend};
        my @extends_new;

        if ($cat_temp->number_of_conformers() > 1) {
            my ($extend_first) = $extend =~ /[Cc]f(\d+)/;
            $extend_first //= 1;
            $extend_first = ($extend_first - 1) * $cf_num + 1;
            @extends_new = map { 'Cf' . $_ } ($extend_first..$extend_first + $cf_num - 1);
        }else {
            @extends_new = ($extend)
        }

        my $new;

        for my $i (0..$#extends_new) {
            my $extend_temp = $extends_new[$i];
            my $new = $newdir;

            my ($extend_partern) = $extend =~ /(\D+)/;
            my ($extend_temp_partern) = $extend_temp =~ /(\D+)/;
            if ($extend_partern ne $extend_temp_partern) {
                $new .= "/$extend_temp";
            }else {
                $new =~ s/$extend/$extend_temp/;
            }

            my $head = $new; $head =~ s/\/Cf\d+//;
            if ( ! exists $status->{$head}) {
                my @pattern = split('/', $head);
                my $state = $pattern[-1];

				$W_Key->{skip_step1} //= $skip_step1;

                if ($state =~ /^ts/i) {
                    $status->{$head} = new Aaron::G09Job_TS(
                        name => $head,
                        catalysis => $new_dir->{$newdir}->{catalysis} ,
                        Gkey => $G_Key,
                        Wkey => $W_Key,
                        template_job => $submission_template,
                        skip_step1 => $W_Key->{skip_step1},
                    );
                }elsif ($state =~ /^min/i) {
                    $status->{$head} = new Aaron::G09Job_MIN(
                        name => $head,
                        catalysis => $new_dir->{$newdir}->{catalysis},
                        Gkey => $G_Key,
                        Wkey => $W_Key,
                        template_job => $submission_template,
                        skip_step1 => $W_Key->{skip_step1},
                    );
                }
            }elsif (%{$status->{$head}->{catalysis}}) {

                $status->{$head}->{catalysis}->conformer_geometry($new_dir->{$newdir}->{catalysis});
            }else {
                $status->{$head}->{catalysis} = $new_dir->{$newdir}->{catalysis};
            }

            if ($extend_temp !~ /Cf/) {
                if (! -d $new && ($status->{$head}->{status} ne 'repeated')) {
                    make_path($new);
                    $status->{$head}->build_com();
                }
            }else {
                if (! $status->{$head}->{conformers}->{$extend_temp} ) {
                    $status->{$head}->new_conformer( conf => $extend_temp );
                }else {
                    $status->{$head}->{conformers}->{$extend_temp}->{catalysis} =
                        $status->{$head}->{catalysis};
                }

                if (! -d $new &&
                    ($status->{$head}->{conformers}->{$extend_temp}->{status} ne 'repeated')) {
                    make_path($new);
                    if ($i == 0 || ($W_Key->{full_conformers})) {
                        my ($num_cf) = $extend_temp =~ /Cf(\d+)/;
                        my $catalysis = $status->{$head}->{conformers}->{$extend_temp}->{catalysis};
                        $catalysis->make_conformer( new_number => $i + 1 );
                        #try to remove clashes in new conformers by wiggling atoms that are too close together
                        my $no_clash = $catalysis->remove_clash();

                        $status->{$head}->{conformers}->{$extend_temp}->build_com();

                        #if the clash was removed, mark this conformer for submission
                        #otherwise, skip it
                        if( $no_clash ) {
                            $status->{$head}->{conformers}->{$extend_temp}->set_status('2submit');
                        } else {
                            $status->{$head}->{conformers}->{$extend_temp}->set_status('skipped');
                            $status->{$head}->{conformers}->{$extend_temp}->set_msg("could not remove clash\n");
                        }
                    }
                }
            }
        }
        $new_dir->{$newdir}->{catalysis}->init_conf_num();
    }
}


sub get_status_run {

    for my $ligand (keys %{ $ligs_subs }) {

        chdir ($ligand);

        foreach my $geometry (sort keys %{ $ligs_subs->{$ligand}->{jobs} }) {
            $ligs_subs->{$ligand}->{jobs}->{$geometry}->check_status_run();
        }

        chdir ($parent);
    }
}


sub analyze_result {

    my $data = '';
    for my $ligand (sort keys %{ $ligs_subs }) {
        $data .= "Available thermochemical data for $ligand (T = $G_Key->{temperature} K):\n";
        my @items = ($ligand, sort keys %{ $ligs_subs->{$ligand}->{substrate} });
        for my $item (@items) {
            $data .= "$item:\n";
            $data .= &_analyze_result($ligs_subs->{$ligand}->{jobs}, $item);
        }
        $data .= "\n";
    }
    print_to_thermo($data);
#    print_message_to_log("Updating thermo file\n");
}


sub _analyze_result {

    my ($jobs, $item) = @_;

    my @geo = grep { $_ =~ /^$item\// } sort keys %{ $jobs };

    #find minimum stereo_geo for each thermo
    my @min = (9999) x 8;
    for my $geo (@geo) {
        for my $i (0..$#min) {
            if (my $conf_hash = $jobs->{$geo}->{conformers}) {
                my @conformers = grep { $conf_hash->{$_}->{thermo}->[$i] }
                                    sort keys %{ $conf_hash };
                my @thermo = sort { $a <=> $b }
                             map {$conf_hash->{$_}->{thermo}->[$i]} @conformers;

                $min[$i] = $thermo[0] < $min[$i] ? $thermo[0] : $min[$i] if @thermo;
            }else {
                $min[$i] = $jobs->{$geo}->{thermo}->[$i] &&
                           ($jobs->{$geo}->{thermo}->[$i] < $min[$i]) ?
                           $jobs->{$geo}->{thermo}->[$i] : $min[$i];
            }
        }
    }

    my $RT = PHYSICAL->{BOLTZMANN} * $G_Key->{temperature};

    my @stereo_geo;

    for my $geo (@{ $W_Key->{selectivity} }) {
        if (my @geo_temp = grep {$_ =~ /\/$geo\//i} @geo) {
            push (@stereo_geo, [@geo_temp]);
        }
    }

    my $no_sele;
    unless (@stereo_geo) {
        @stereo_geo = ([@geo]);
        $no_sele = 1;
    }

    my $thermo = {};
	my $abs_thermo = {};

    for my $n (0..$#stereo_geo) {
        my $key = $no_sele ? 'NONE' : $W_Key->{selectivity}->[$n];
        $thermo->{$key} = {};
        $thermo->{$key}->{sum} = [];
        $thermo->{$key}->{geos} = {};

        $abs_thermo->{$key} = {};
        $abs_thermo->{$key}->{geos} = {};

        my @geos = @{ $stereo_geo[$n] };

        for my $geo (@geos) {
            if ($jobs->{$geo}->{conformers}) {
                $thermo->{$key}->{geos}->{$geo}->{conformers} = {};
                $abs_thermo->{$key}->{geos}->{$geo}->{conformers} = {};
                my $thermo_cf_exp = [];

                for my $cf (sort keys %{ $jobs->{$geo}->{conformers} }) {
                    my $job = $jobs->{$geo}->{conformers}->{$cf};

                    ! @{$job->{thermo}} && do { next; };

                    my @thermo_rel;
                    if (! $W_Key->{multistep}) {
                        @thermo_rel = map { ($job->{thermo}->[$_] - $min[$_]) * $hart_to_kcal }
                                          (0..$#{ $job->{thermo} });
                    }

                    for my $i (0..$#thermo_rel) {
                        $thermo_cf_exp->[$i] += exp(-$thermo_rel[$i]/$RT);
                    }

                    $thermo->{$key}->{geos}->{$geo}->{conformers}->{$cf} = [@thermo_rel] if @thermo_rel;
                    $abs_thermo->{$key}->{geos}->{$geo}->{conformers}->{$cf} = [@{ $job->{thermo} }];
                    $abs_thermo->{$key}->{found} = 1;
                }
                my @thermo_cf = map { -1*$RT*log($_) } @$thermo_cf_exp if @$thermo_cf_exp;
                $thermo->{$key}->{geos}->{$geo}->{thermo} = [@thermo_cf] if @$thermo_cf_exp;
            }else {
                for my $i (0..$#{ $jobs->{$geo}->{thermo} }) {
                    if (! $W_Key->{multistep}) {
                        $thermo->{$key}->{geos}->{$geo}->{thermo}->[$i] =
                            ($jobs->{$geo}->{thermo}->[$i] - $min[$i]) * $hart_to_kcal;
                    }
                }
                $abs_thermo->{$key}->{geos}->{$geo}->{thermo} =  [@{ $jobs->{$geo}->{thermo} }];
                if (@{ $jobs->{$geo}->{thermo} } ) {
                    $thermo->{$key}->{found} = 1;
                    $abs_thermo->{$key}->{found} = 1;
                }
            }

            for my $i (0.. $#{ $thermo->{$key}->{geos}->{$geo}->{thermo} }) {
                $thermo->{$key}->{sum}->[$i] +=
                    exp(-$thermo->{$key}->{geos}->{$geo}->{thermo}->[$i]/$RT);
            }
        }
    }

    my $data = "Relative thermochemistry (kcal/mol)\n";
    $data .= print_ee($thermo);

    $data .= "Absolute thermochemistry (hartees)\n";
    $data .= print_ee($abs_thermo, 1);

    return $data;
}


sub get_flags {
    my ($file) = @_;

    my @flags;
    my @coords = grab_coords($file);

    for my $atom (0..$#coords) {
        if ($coords[$atom][1] == 0) {
            push (@flags, $atom);
        }
    }

    return [@flags];
}

1;
