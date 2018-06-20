#!/usr/bin/perl -w

=head1 SYNOPSIS
Provides utilities for common interactions with AaronTools
=cut

use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

package _utils;
use strict;
use Data::Dumper;

sub get_geom {

=item get_geom($file)
Gets geometry object from file
=cut

    use AaronTools::Geometry;

    my $file = shift;
    my $geom = new AaronTools::Geometry();
    $geom->{name} = $file;
    $geom->read_geometry($file);
    unless ( @{ $geom->elements() } ) {
        print {*STDERR} ("\nCouldn't read geometry: $file\n\n");
        return 0;
    }
    return $geom;
}

sub get_cat {

=item get_cat($file, \%substituents)
Gets catalysis object from file. Optional substituent information can be provided.

\%substituents = {'ligand'=>{ atom=>label, ... }, 'substrate'=>{ atom=>label, ... }}
=cut

    use AaronTools::Catalysis;

    my $file = shift;

    my %params;
    $params{substituents} = shift;
    ( $params{name} ) = ( $file =~ /(.*)\..*?$/ );

    my $relnum = shift;

    # creat Catalysis object
    my $cat = new AaronTools::Catalysis( 'name' => $params{name} );

    # if substituent info provided, make atom numbering relative, if necessary
    if ( $params{substituents} && !($relnum) ) {
        # this gets the number to subtract for relative numbering
        my $ligstart = get_ligstart($cat);
        my $substart = get_substart($cat);

        # update info hash to reflect numbering change
        foreach my $sub ( keys %{ $params{substituents}{ligand} } ) {
            $params{substituents}{ligand}{ $sub - $ligstart } =
              delete $params{substituents}{ligand}{$sub};
        }
        foreach my $sub ( keys %{ $params{substituents}{substrate} } ) {
            $params{substituents}{substrate}{ $sub - $substart } =
              delete $params{substituents}{substrate}{$sub};
        }

    }
    # update catalysis object with substituent info, if provided
    $cat = new AaronTools::Catalysis(%params) if ( $params{substituents} );

    # check for success
    unless ( @{ $cat->{elements} } ) {
        print {*STDERR} ("\nCouldn't read catalyst geometry: $file\n\n");
        return 0;
    }
    return $cat;
}

sub get_lig {

=item get_lig($file)
Reads in ligand object from .xyz file or by built-in name
=cut

    use AaronTools::Catalysis;

    my $file = shift;
    my $lig;
    if ( $file =~ /.*\.xyz$/ ) {
        $lig = new AaronTools::Ligand( name => ( $file =~ /(.*)\..*?$/ ) );
    } else {
        $lig = new AaronTools::Ligand( name => $file );
    }
    unless ( @{ $lig->{elements} } ) {
        print {*STDERR} ("\nCouldn't read ligand geometry: $file\n\n");
        return 0;
    }
    return $lig;
}

sub get_outfile {

=item get_outfile($filebase, $path, \@appends, $sep)
Generates an outfile name for printXYZ() methods

outfile = path/filebase_appends.xyz (sep defaults to _)
=cut

    # prints to STDOUT if $path == ''
    # or saves to infile_append1_append2_etc.xyz
    # $path= '-', defaults to cwd
    my $filebase = shift;
    my $path     = shift;
    my $appends  = shift(@_) // [];
    my $sep      = '_';

    my $outfile = '';
    if ( $path ne '-' ) {
        # strip just file name (no path or file extension)
        $outfile = $filebase;
        $outfile =~ s/(.*\/)?(.*)\..*?$/$2/;
        if ( $path ne '' ) {
            # if no directory specified, write to cwd
            # make sure we don't have double path seperators!
            if ( $path =~ /.*\/$/ ) {
                $outfile = $path . $outfile;
            } else {
                $outfile = $path . '/' . $outfile;
            }
        }
        foreach my $append (@$appends) {
            $outfile .= $sep . $append;
        }
        $outfile .= '.xyz';
    }
    return $outfile;
}

sub strip_dir {

=item strip_dir($fname)
Removes the directory path and returns only the file name
=cut

    my $fname = shift;
    $fname =~ s/.*\/(.*)/$1/;
    return $fname;
}

sub get_ligstart {

=item get_ligstart($catalysis)
Returns a value to be subtracted from an atom index to switch from absolute to relative indexing
=cut

    my $cat = shift;
    return ( sort { $a <=> $b } @{ $cat->{ligand_atoms} } )[0];
}

sub get_substart {

=item get_substart($catalysis)
Returns a value to be subtracted from an atom index to switch from absolute to relative indexing
=cut

    my $cat = shift;
    return ( sort { $a <=> $b } @{ $cat->{substrate_atoms} } )[0];
}

# Used to determine AaronTools atom numbering from command-line provided numbering
sub get_sub_numbers {
	my $geom      = shift;
	my $component = shift;

	my %convert;
	# determine requested target atom, in aaron's number scheme
	foreach my $aaron_num ( keys %{ $geom->{$component}{substituents} } ) {
		my $target;
		if ( $geom->{$component}->{substituents}->{$aaron_num}->{sub} ) {
			# only add to screen_sub convert if it was a substituent requested
			$target =
				$geom->{$component}->{substituents}->{$aaron_num}->{sub};
		} else {
			# but not ones that were simply auto-detected
			next;
		}

		$convert{$target} = $aaron_num;
	}
	return %convert;
}

1;
