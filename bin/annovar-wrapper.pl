#!/usr/bin/env perl 

package Main;

use lib "/home/guests/jir2004/perlmodule/Annovar-Wrapper/lib";

use Moose;
#use Carp::Always;

extends 'Annovar::Wrapper';

Annovar::Wrapper->new_with_options->run;

1;
