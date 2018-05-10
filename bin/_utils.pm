#!/usr/bin/perl -w

use strict;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

use AaronTools::Geometry;

package _utils;


sub get_outfile {

    # prints to STDOUT if $path == ''
    # or saves to infile_append1_append2_etc.xyz
    # $path= '-', defaults to cwd
    my $filebase = shift;
    my $path     = shift;
    my $appends  = shift(@_) // [];
	my $sep = '_';

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

sub get_geom {
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

sub strip_dir {
	my $fname = shift;
	$fname =~ s/.*\/(.*)/$1/;
	return $fname;
}

1;
