package AaronTools::Molecules;

use strict;
use Exporter qw(import);

our @EXPORT = qw(built_in coords elements flags);

my $molecules = { methane => 
                    { coords=>[[0.92646611, 1.86569721, -0.87365150],
                               [1.28312053, 0.85688721, -0.87365150],
                               [1.28313895, 2.37009540, 0.00000000],
                               [1.28313895, 2.37009540, -1.74730301],
                               [-0.14353389, 1.86571040, -0.87365150]],
                      elements=>['C', 'H', 'H', 'H', 'H'],
                      flags => [0, 0, 0, 0, 0],
                    },
                   benzene =>
                    { coords=>[[-0.70156658, -1.21629378, -0.00045000], 
                               [ 0.69359342, -1.21629378, -0.00045000], 
                               [ 1.39113142, -0.00854278, -0.00045000], 
                               [ 0.69347742,  1.19996622, -0.00164900], 
                               [-0.70134758,  1.19988822, -0.00212800], 
                               [-1.39894858, -0.00831778, -0.00113200], 
                               [-1.25132558, -2.16861078,  0.00000000], 
                               [ 1.24310142, -2.16880678,  0.00086500], 
                               [ 2.49081142, -0.00846278,  0.00018400], 
                               [ 1.24367742,  2.15210922, -0.00170800], 
                               [-1.25146958,  2.15216922, -0.00308100], 
                               [-2.49855258, -0.00813478, -0.00131200]],
                      elements=>[qw(C C C C C C H H H H H H)],
                      flags=>[0,0,0,0,0,0,0,0,0,0,0,0],
                    } 
   
                 };

sub built_in {
    my $name = shift;
    if (exists $molecules->{$name}) {
        return 1;
    }else {
        return 0;
    }
}


sub coords {
    my $name = shift;
    if (built_in($name)) {
        return $molecules->{$name}->{coords};
    }else {
        return [];
    }
}


sub elements {
    my $name = shift;
    if (built_in($name)) {
        return $molecules->{$name}->{elements};
    }else {
        return [];
    }
}


sub flags {
    my $name = shift;
    if (built_in($name)) {
        return $molecules->{$name}->{flags};
    }else {
        return [];
    }
}

1;
