#!/usr/bin/env perl 

package Main;

use Moose;
#use Carp::Always;

extends 'Annovar::Wrapper';

Annovar::Wrapper->new_with_options->run;

1;
