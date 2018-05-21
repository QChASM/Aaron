package Math::Vector::Real;

our $VERSION = '0.10';

use strict;
use warnings;
use Carp;
use POSIX ();

use Exporter qw(import);
our @EXPORT = qw(V);

local ($@, $!, $SIG{__DIE__});
eval { require Math::Vector::Real::XS };

our %op = (add => '+',
	   neg => 'neg',
	   sub => '-',
	   mul => '*',
	   div => '/',
	   cross => 'x',
	   add_me => '+=',
	   sub_me => '-=',
	   mul_me => '*=',
	   div_me => '/=',
	   abs => 'abs',
	   atan2 => 'atan2',
	   equal => '==',
	   nequal => '!=',
	   clone => '=',
	   as_string => '""');

our %ol;
$ol{$op{$_}} = \&{${Math::Vector::Real::}{$_}} for keys %op;

require overload;
overload->import(%ol);

sub V { bless [@_] }

sub new {
    my $class = shift;
    bless [@_], $class
}

sub new_ref {
    my $class = shift;
    bless [@{shift()}], $class;
}

sub zero {
    my ($class, $dim) = @_;
    $dim >= 0 or croak "negative dimension";
    bless [(0) x $dim], $class
}

sub is_zero {
    $_ and return 0 for @$_[0];
    return 1
}

sub cube {
    my ($class, $dim, $size) = @_;
    bless [($size) x $dim], $class;
}

sub axis_versor {
    my ($class, $dim, $ix);
    if (ref $_[0]) {
        my ($self, $ix) = @_;
        $class = ref $self;
        $dim = @$self;
    }
    else {
        ($class, $dim, $ix) = @_;
        $dim >= 0 or croak "negative dimension";
    }
    ($ix >= 0 and $ix < $dim) or croak "axis index out of range";
    my $self = [(0) x $dim];
    $self->[$ix] = 1;
    bless $self, $class
}

sub _caller_op {
    my $level = (shift||1) + 1;
    my $sub = (caller $level)[3];
    $sub =~ s/.*:://;
    my $op = $op{$sub};
    (defined $op ? $op : $sub);
}

sub _check_dim {
    local ($@, $SIG{__DIE__});
    eval { @{$_[0]} == @{$_[1]} } and return;
    my $op = _caller_op(1);
    my $loc = ($_[2] ? 'first' : 'second');
    UNIVERSAL::isa($_[1], 'ARRAY') or croak "$loc argument to vector operator '$op' is not a vector";
    croak "vector dimensions do not match";
}

sub clone { bless [@{$_[0]}] }

sub set {
    &_check_dim;
    my ($v0, $v1) = @_;
    $v0->[$_] = $v1->[$_] for 0..$#$v1;
}

sub add {
    &_check_dim;
    my ($v0, $v1) = @_;
    bless [map $v0->[$_] + $v1->[$_], 0..$#$v0]
}

sub add_me {
    &_check_dim;
    my ($v0, $v1) = @_;
    $v0->[$_] += $v1->[$_] for 0..$#$v0;
    $v0;
}

sub neg { bless [map -$_, @{$_[0]}] }

sub sub {
    &_check_dim;
    my ($v0, $v1) = ($_[2] ? @_[1, 0] : @_);
    bless [map $v0->[$_] - $v1->[$_], 0..$#$v0]
}

sub sub_me {
    &_check_dim;
    my ($v0, $v1) = @_;
    $v0->[$_] -= $v1->[$_] for 0..$#$v0;
    $v0;
}

sub mul {
    if (ref $_[1]) {
	&_check_dim;
	my ($v0, $v1) = @_;
	my $acu = 0;
	$acu += $v0->[$_] * $v1->[$_] for 0..$#$v0;
	$acu;
    }
    else {
	my ($v, $s) = @_;
	bless [map $s * $_, @$v];
    }
}

sub mul_me {
    ref $_[1] and croak "can not multiply by a vector in place as the result is not a vector";
    my ($v, $s) = @_;
    $_ *= $s for @$v;
    $v
}

sub div {
    croak "can't use vector as dividend"
	if ($_[2] or ref $_[1]);
    my ($v, $div) = @_;
    $div == 0 and croak "illegal division by zero";
    my $i = 1 / $div;
    bless [map $i * $_, @$v]
}

sub div_me {
    croak "can't use vector as dividend" if ref $_[1];
    my $v = shift;
    my $i = 1 / shift;
    $_ *= $i for @$v;
    $v;
}

sub equal {
    &_check_dim;
    my ($v0, $v1) = @_;
    $v0->[$_] == $v1->[$_] || return 0 for 0..$#$v0;
    1;
}

sub nequal {
    &_check_dim;
    my ($v0, $v1) = @_;
    $v0->[$_] == $v1->[$_] || return 1 for 0..$#$v0;
    0;
}

sub cross {
    &_check_dim;
    my ($v0, $v1) = ($_[2] ? @_[1, 0] : @_);
    my $dim = @$v0;
    if ($dim == 3) {
	return bless [$v0->[1] * $v1->[2] - $v0->[2] * $v1->[1],
		      $v0->[2] * $v1->[0] - $v0->[0] * $v1->[2],
		      $v0->[0] * $v1->[1] - $v0->[1] * $v1->[0]]
    }
    if ($dim == 7) {
	croak "cross product for dimension 7 not implemented yet, patches welcome!";
    }
    else {
	croak "cross product not defined for dimension $dim"
    }
}

sub as_string { "{" . join(", ", @{$_[0]}). "}" }

sub abs {
    my $acu = 0;
    $acu += $_ * $_ for @{$_[0]};
    sqrt $acu;
}

*norm = \&abs;

sub abs2 {
    my $acu = 0;
    $acu += $_ * $_ for @{$_[0]};
    $acu;
}

*norm2 = \&abs2;

sub dist {
    &_check_dim;
    my ($v0, $v1) = @_;
    my $d2 = 0;
    for (0..$#$v0) {
	my $d = $v0->[$_] - $v1->[$_];
	$d2 += $d * $d;
    }
    sqrt($d2);
}

sub dist2 {
    &_check_dim;
    my ($v0, $v1) = @_;
    my $d2 = 0;
    for (0..$#$v0) {
	my $d = $v0->[$_] - $v1->[$_];
	$d2 += $d * $d;
    }
    $d2;
}

sub manhattan_norm {
    my $n = 0;
    $n += CORE::abs($_) for @{$_[0]};
    return $n;
}

sub manhattan_dist {
    &_check_dim;
    my ($v0, $v1) = @_;
    my $d = 0;
    $d += CORE::abs($v0->[$_] - $v1->[$_]) for 0..$#$v0;
    return $d;
}

sub _upgrade {
    my $dim;
    map {
	my $d = eval { @{$_} };
	defined $d or croak "argument is not a vector or array";
	if (defined $dim) {
	    $d == $dim or croak "dimensions do not match";
	}
	else {
	    $dim = $d;
	}
	UNIVERSAL::isa($_, __PACKAGE__) ? $_ : clone($_);
    } @_;
}

sub atan2 {
    my ($v0, $v1) = @_;
    my $a0 = &abs($v0);
    return 0 unless $a0;
    my $u0 = $v0 / $a0;
    my $p = $v1 * $u0;
    CORE::atan2(&abs($v1 - $p * $u0), $p);
}

sub versor {
    my $self = shift;
    my $f = 0;
    $f += $_ * $_ for @$self;
    $f == 0 and croak "Illegal division by zero";
    $f = 1/sqrt $f;
    bless [map $f * $_, @$self]
}

sub wrap {
    my ($self, $v) = @_;
    &_check_dim;

    bless [map  { my $s = $self->[$_];
		  my $c = $v->[$_];
		  $c - $s * POSIX::floor($c/$s) } (0..$#$self)];
}

sub max_component {
    my $max = 0;
    for (@{shift()}) {
	my $abs = CORE::abs($_);
	$abs > $max and $max = $abs;
    }
    $max
}

sub min_component {
    my $self = shift; 
    my $min = CORE::abs($self->[0]);
    for (@$self) {
	my $abs = CORE::abs($_);
	$abs < $min and $min = $abs;
    }
    $min
}

*max = \&max_component;
*min = \&min_component;

sub box {
    shift;
    return unless @_;
    my $min = clone(shift);
    my $max = clone($min);
    my $dim = $#$min;
    for (@_) {
        for my $ix (0..$dim) {
            my $c = $_->[$ix];
            if ($max->[$ix] < $c) {
                $max->[$ix] = $c;
            }
            elsif ($min->[$ix] > $c) {
                $min->[$ix] = $c
            }
        }
    }
    ($min, $max);
}

sub max_component_index {
    my $self = shift;
    return unless @$self;
    my $max = 0;
    my $max_ix = 0;
    for my $ix (0..$#$self) {
        my $c = CORE::abs($self->[$ix]);
        if ($c > $max) {
            $max = $c;
            $max_ix = $ix;
        }
    }
    $max_ix;
}

sub min_component_index {
    my $self = shift;
    return unless @$self;
    my $min = CORE::abs($self->[0]);
    my $min_ix = 0;
    for my $ix (1..$#$self) {
        my $c = CORE::abs($self->[$ix]);
        if ($c < $min) {
            $min = $c;
            $min_ix = $ix
        }
    }
    $min_ix;
}

sub decompose {
    my ($u, $v) = @_;
    my $p = $u * ($u * $v)/abs2($u);
    my $n = $v - $p;
    wantarray ? ($p, $n) : $n;
}

sub canonical_base {
    my ($class, $dim) = @_;
    my @base = map { bless [(0) x $dim], $class } 1..$dim;
    $base[$_][$_] = 1 for 0..$#base;
    return @base;
}

sub rotation_base_3d {
    my $v = shift;
    @$v == 3 or croak "rotation_base_3d requires a vector with three dimensions";
    $v = $v->versor;
    my $n = [0, 0, 0];
    for (0..2) {
        if (CORE::abs($v->[$_]) > 0.57) {
            $n->[($_ + 1) % 3] = 1;
            $n = $v->decompose($n)->versor;
            return ($v, $n, $v x $n);
        }
    }
    die "internal error, all the components where smaller than 0.57!";
}

sub rotate_3d {
    my $v = shift;
    my $angle = shift;
    my $c = cos($angle); my $s = sin($angle);
    my ($i, $j, $k) = $v->rotation_base_3d;
    my $rj = $c * $j + $s * $k;
    my $rk = $c * $k - $s * $j;
    if (wantarray) {
        return map { ($_ * $i) * $i + ($_ * $j) * $rj + ($_ * $k) * $rk } @_;
    }
    else {
        my $a = shift;
        return (($a * $i) * $i + ($a * $j) * $rj + ($a * $k) * $rk);
    }
}

sub normal_base { __PACKAGE__->complementary_base(@_) }

sub complementary_base {
    shift;
    @_ or croak "complementaty_base requires at least one argument in order to determine the dimension";
    my $dim = @{$_[0]};
    if ($dim == 2 and @_ == 1) {
        my $u = versor($_[0]);
        @$u = ($u->[1], -$u->[0]);
        return $u;
    }

    my @v = map clone($_), @_;
    my @base = Math::Vector::Real->canonical_base($dim);
    for my $i (0..$#v) {
        my $u = versor($v[$i]);
        $_ = decompose($u, $_) for @v[$i+1 .. $#v];
        $_ = decompose($u, $_) for @base;
    }

    my $last = $#base - @v;
    return if $last < 0;
    for my $i (0 .. $last) {
        my $max = abs2($base[$i]);
        if ($max < 0.3) {
            for my $j ($i+1 .. $#base) {
                my $d2 = abs2($base[$j]);
                if ($d2 > $max) {
                    @base[$i, $j] = @base[$j, $i];
                    last unless $d2 < 0.3;
                    $max = $d2;
                }
            }
        }
        my $versor = $base[$i] = versor($base[$i]);
        $_ = decompose($versor, $_) for @base[$i+1..$#base];
    }
    wantarray ? @base[0..$last] : $base[0];
}

sub select_in_ball {
    my $v = shift;
    my $r = shift;
    my $r2 = $r * $r;
    grep $v->dist2($_) <= $r2, @_;
}

sub select_in_ball_ref2bitmap {
    my $v = shift;
    my $r = shift;
    my $p = shift;
    my $r2 = $r * $r;
    my $bm = "\0" x int((@$p + 7) / 8);
    for my $ix (0..$#$p) {
        vec($bm, $ix, 1) = 1 if $v->dist2($p->[$ix]) <= $r2;
    }
    return $bm;
}

1;
__END__

=head1 NAME

Math::Vector::Real - Real vector arithmetic in Perl

=head1 SYNOPSIS

  use Math::Vector::Real;

  my $v = V(1.1, 2.0, 3.1, -4.0, -12.0);
  my $u = V(2.0, 0.0, 0.0,  1.0,   0.3);

  printf "abs(%s) = %d\n", $v, abs($b);
  my $dot = $u * $v;
  my $sub = $u - $v;
  # etc...

=head1 DESCRIPTION

A simple pure perl module to manipulate vectors of any dimension.

The function C<V>, always exported by the module, allows one to create
new vectors:

  my $v = V(0, 1, 3, -1);

Vectors are represented as blessed array references. It is allowed to
manipulate the arrays directly as far as only real numbers are
inserted (well, actually, integers are also allowed because from a
mathematical point of view, integers are a subset of the real
numbers).

Example:

  my $v = V(0.0, 1.0);

  # extending the 2D vector to 3D:
  push @$v, 0.0;

  # setting some component value:
  $v->[0] = 23;

Vectors can be used in mathematical expressions:

  my $u = V(3, 3, 0);
  $p = $u * $v;       # dot product
  $f = 1.4 * $u + $v; # scalar product and vector addition
  $c = $u x $v;       # cross product, only defined for 3D vectors
  # etc.

The currently supported operations are:

  + * /
  - (both unary and binary)
  x (cross product for 3D vectors)
  += -= *= /= x=
  == !=
  "" (stringfication)
  abs (returns the norm)
  atan2 (returns the angle between two vectors)

That, AFAIK, are all the operations that can be applied to vectors.

When an array reference is used in an operation involving a vector, it
is automatically upgraded to a vector. For instance:

  my $v = V(1, 2);
  $v += [0, 2];

=head2 Extra methods

Besides the common mathematical operations described above, the
following methods are available from the package.

Note that all these methods are non destructive returning new objects
with the result.

=over 4

=item $v = Math::Vector::Real->new(@components)

Equivalent to C<V(@components)>.

=item $zero = Math::Vector::Real->zero($dim)

Returns the zero vector of the given dimension.

=item $v = Math::Vector::Real->cube($dim, $size)

Returns a vector of the given dimension with all its components set to
C<$size>.

=item $u = Math::Vector::Real->axis_versor($dim, $ix)

Returns a unitary vector of the given dimension parallel to the axis
with index C<$ix> (0-based).

For instance:

  Math::Vector::Real->axis_versor(5, 3); # V(0, 0, 0, 1, 0)
  Math::Vector::Real->axis_versor(2, 0); # V(1, 0)

=item @b = Math::Vector::Real->canonical_base($dim)

Returns the canonical base for the vector space of the given
dimension.

=item $u = $v->versor

Returns the versor for the given vector.

It is equivalent to:

  $u = $v / abs($v);

=item $wrapped = $w->wrap($v)

Returns the result of wrapping the given vector in the box
(hyper-cube) defined by C<$w>.

Long description:

Given the vector C<W> and the canonical base C<U1, U2, ...Un> such
that C<W = w1*U1 + w2*U2 +...+ wn*Un>. For every component C<wi> we
can consider the infinite set of affine hyperplanes perpendicular to
C<Ui> such that they contain the point C<j * wi * Ui> being C<j> an
integer number.

The combination of all the hyperplanes defined by every component
define a grid that divides the space into an infinite set of affine
hypercubes. Every hypercube can be identified by its lower corner
indexes C<j1, j2, ..., jN> or its lower corner point C<j1*w1*U1 +
j2*w2*U2 +...+ jn*wn*Un>.

Given the vector C<V>, wrapping it by C<W> is equivalent to finding
where it lays relative to the lower corner point of the hypercube
inside the grid containing it:

  Wrapped = V - (j1*w1*U1 + j2*w2*U2 +...+ jn*wn*Un)

  such that ji*wi <= vi <  (ji+1)*wi

=item $max = $v->max_component

Returns the maximum of the absolute values of the vector components.

=item $min = $v->min_component

Returns the minimum of the absolute values of the vector components.

=item $d2 = $b->norm2

Returns the norm of the vector squared.

=item $d = $v->dist($u)

Returns the distance between the two vectors.

=item $d = $v->dist2($u)

Returns the distance between the two vectors squared.

=item ($bottom, $top) = Math::Vector::Real->box($v0, $v1, $v2, ...)

Returns the two corners of a hyper-box containing all the given
vectors.

=item $v->set($u)

Equivalent to C<$v = $u> but without allocating a new object.

Note that this method is destructive.

=item $d = $v->max_component_index

Return the index of the vector component with the maximum size.

=item ($p, $n) = $v->decompose($u)

Decompose the given vector C<$u> in two vectors: one parallel to C<$v>
and another normal.

In scalar context returns the normal vector.

=item @b = Math::Vector::Real->complementary_base(@v)

Returns a base for the subspace complementary to the one defined by
the base @v.

The vectors on @v must be linearly independent. Otherwise a division
by zero error may pop up or probably due to rounding errors, just a
wrong result may be generated.

=item @b = $v->normal_base

Returns a set of vectors forming an ortonormal base for the hyperplane
normal to $v.

In scalar context returns just some unitary vector normal to $v.

Note that this two expressions are equivalent:

  @b = $v->normal_base;
  @b = Math::Vector::Real->complementary_base($v);

=item ($i, $j, $k) = $v->rotation_base_3d

Given a 3D vector, returns a list of 3 vectors forming an orthonormal
base where $i has the same direction as the given vector C<$v> and
C<$k = $i x $j>.

=item @r = $v->rotate_3d($angle, @s)

Returns the vectors C<@u> rotated around the vector C<$v> an
angle C<$angle> in radians in anticlockwise direction.

See L<http://en.wikipedia.org/wiki/Rotation_operator_(vector_space)>.

=item @s = $center->select_in_ball($radius, $v1, $v2, $v3, ...)

Selects from the list of given vectors those that lay inside the
n-ball determined by the given radius and center (C<$radius> and
C<$center> respectively).

=back

=head2 Zero vector handling

Passing the zero vector to some methods (i.e. C<versor>, C<decompose>,
C<normal_base>, etc.) is not acceptable. In those cases, the module
will croak with an "Illegal division by zero" error.

C<atan2> is an exceptional case that will return 0 when any of its
arguments is the zero vector (for consistency with the C<atan2> builtin
operating over real numbers).

In any case note that, in practice, rounding errors frequently cause
the check for the zero vector to fail resulting in numerical
instabilities.

The correct way to handle this problem is to introduce in your code
checks of this kind:

  if ($v->norm2 < $epsilon2) {
    croak "$v is too small";
  }

Or even better, reorder the operations to minimize the chance of
instabilities if the algorithm allows it.

=head1 SEE ALSO

L<Math::Vector::Real::Random> extends this module with random vector
generation methods.

L<Math::GSL::Vector>, L<PDL>.

There are other vector manipulation packages in CPAN (L<Math::Vec>,
L<Math::VectorReal>, L<Math::Vector>), but they can only handle 3
dimensional vectors.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009-2012 by Salvador FandiE<ntilde>o
(sfandino@yahoo.com)

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=cut
