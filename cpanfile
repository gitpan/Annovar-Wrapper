requires 'Cwd';
requires 'Data::Dumper';
requires 'DateTime';
requires 'File::Basename';
requires 'File::Find::Rule';
requires 'File::FindLib';
requires 'File::Path';
requires 'IO::Uncompress::Gunzip';
requires 'Moose';
requires 'MooseX::Getopt';
requires 'perl', '5.006';

on build => sub {
    requires 'Test::More';
};
