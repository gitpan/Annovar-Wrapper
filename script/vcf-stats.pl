#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
# Completely stole this file - just wanted to add a method to validate each line as it passes
#

use strict;
use warnings;
use Carp;
use VcfStats;
use Data::Dumper;

my $opts = parse_params();
vcf_stats($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg )
    {
        croak @msg;
    }
    die
    "Usage: vcf-stats [OPTIONS] file.vcf.gz\n",
    "Options:\n",
    "   -d, --dump <file>                           Take an existing dump file and recreate the files (works with -p)\n",
    "   -f, --filters <filter1,filter2>             List of filters such as column/field (any value), column/field=bin:max (cluster in bins),column/field=value (exact value)\n",
    "   -p, --prefix <dir/string>                   Prefix of output files. If slashes are present, directories will be created.\n",
    "   -s, --samples <list>                        Process only the listed samples, - for none. Excluding unwanted samples may increase performance considerably.\n",
    "   -h, -?, --help                              This help message.\n",
    "\n",
    "Examples:\n",
    "   # Calculate stats separately for the filter field, quality and non-indels\n",
    "   vcf-stats file.vcf.gz -f FILTER,QUAL=10:200,INFO/INDEL=False -p out/\n",
    "\n",
    "   # Calculate stats for all samples\n",
    "   vcf-stats file.vcf.gz -f FORMAT/DP=10:200 -p out/\n",
    "\n",
    "   # Calculate stats only for the sample NA00001\n",
    "   vcf-stats file.vcf.gz -f SAMPLE/NA00001/DP=1:200 -p out/\n",
    "\n",
    "   vcf-stats file.vcf.gz > perl.dump\n",
    "\n";
}


sub parse_params
{
    my $opts = { filters=>{}, filter_param=>'' };
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-d' || $arg eq '--dump'  ) { $$opts{dump}=shift(@ARGV); next; }
        if ( $arg eq '-f' || $arg eq '--filters'  ) { $$opts{filter_param}=shift(@ARGV); next; }
        if ( $arg eq '-p' || $arg eq '--prefix'  ) { $$opts{prefix}=shift(@ARGV); next; }
        if ( $arg eq '-s' || $arg eq '--samples'  ) 
        { 
            my $samples = shift(@ARGV);
            $$opts{samples} = [ split(/,/,$samples) ];
            next;
        }
        if ( -e $arg ) { $$opts{file} = $arg; next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter or nonexistent file: \"$arg\". Run -h for help.\n");
    }
    if ( exists($$opts{dump}) && !exists($$opts{prefix}) ) { error("Expected -p option with -d.\n"); }
    return $opts;
}


sub init_filters
{
    my ($opts,$vcf) = @_;

    for my $filter (split(/,/,$$opts{filter_param}))
    {
        my ($key,$value) = split(/=/,$filter);

        my $rec = { value=>$value, exact=>0, any=>0, bin=>0, is_flag=>0 };
        if ( $key=~m{^INFO/} )
        {
            my $tag = $';
            $$rec{tag} = $tag;
            if ( exists($$vcf{header}{'INFO'}) && exists($$vcf{header}{'INFO'}{$tag}) && $$vcf{header}{'INFO'}{$tag}{Type} eq 'Flag' )
            {
                $$rec{is_flag} = 1;
                $$rec{value} = $value eq 'False' ? 0 : 1;
                $key = "INFO/$tag=". ($$rec{value} ? 'True':'False');
            }
        }
        elsif ( $key eq 'INFO' )
        {
            # All INFO flags should be counted
            for my $tag (keys %{$$vcf{header}{'INFO'}})
            {
                if ( $$vcf{header}{'INFO'}{$tag}{Type} ne 'Flag' ) { next; }
                $$opts{filters}{"INFO/$tag=True"} = { %$rec, is_flag=>1, value=>1, tag=>$tag };
            }
            next;
        }

        if ( ! defined $value )
        {
            $$rec{any} = 1;
        }
        elsif ( $value=~/^(.+):(.+)$/ ) 
        {
            $$rec{bin}      = 1;
            $$rec{bin_size} = $1;
            $$rec{max}      = $2;
        }
        else
        {
            $$rec{exact} = 1;
        }
        $$opts{filters}{$key} = $rec;
    }
}


sub vcf_stats
{
    my ($opts) = @_;

    my $warnref = {
        Header => 0,
        EmptyColname => 0,
        MalformedPos => 0,
        MalformedChr => 0,
        Sorting => 0,
        RefBase => 0,
        Duplicate => 0,
        ALT_field => 0,
        QUAL_field => 0,
        FILTER_field => 0,
        INFO_field => 0,
        Genotypes => 0,
        Total => 0,
        TotalLinesError => 0,
        TotalLines => 0,
        line => 0,
    };

    if ( exists($$opts{dump}) )
    {
        # Use existing dump to recreate the files
        my $vcf = VcfStats->new(file=>'/dev/null');
        $$vcf{stats} = do $$opts{dump};
        $vcf->save_stats($$opts{prefix});
        return;
    }

    # Open the VCF file
    my $vcf = $$opts{file} ? VcfStats->new(file=>$$opts{file}) : VcfStats->new(fh=>\*STDIN);
    init_filters($opts,$vcf);

    $warnref = begin_validate($vcf, $warnref);

    my (@samples) = $vcf->get_samples();

    # Include only requested samples
    if ( exists $$opts{samples} )
    {   
        my @include = ();
        if ( scalar @{$$opts{samples}}>1 or $$opts{samples}[0] ne '-' ) 
        { 
            for my $sample (@{$$opts{samples}}) { push @include,$sample; }
        }
        $vcf->set_samples(include=>\@include); 
    }

    while (my $line=$vcf->next_line())
    {
        next unless $line;
        $warnref->{TotalLines} +=1;
        my $aref = $vcf->next_data_array($line);
        $warnref = validate($vcf, $aref, $warnref);
        if($warnref->{line}){
            $warnref->{Total} += $warnref->{line};
            $warnref->{line} = 0;
            $warnref->{TotalLinesError} += 1;
            next;
        }
        my $href=$vcf->next_data_hash($line);
        $vcf->collect_stats($href,$$opts{filters});
    }

    $vcf->save_stats($$opts{prefix});
#    parse_dump($vcf);
}

sub parse_dump{
    my($self) = @_;
    print Dumper($self->dump);
}

sub begin_validate{
    my($self, $warn) = @_;

    $self->parse_header();
    $self->validate_header();

    if ( !exists($$self{header}) ) {
        $self->warn(qq[The header not present.\n]);
        $warn->{Header} += 1;
    }
    elsif ( !exists($$self{header}{fileformat}) ){
        $warn->{Header} += 1;
        $self->warn(qq[The "fileformat" field not present in the header, assuming VCFv$$self{version}\n]);
    }
    elsif ( $$self{header_lines}[0]{key} ne 'fileformat' ) {
        $warn->{Header} += 1;
        $self->warn(qq[The "fileformat" not the first line in the header\n]);
    }
    if ( !exists($$self{columns}) ) { 
        $self->warn("No column descriptions found.\n"); 
        $warn->{Header} += 1;
    }

    return $warn;
}

sub validate{
    my ($self, $line, $warn) = @_;

    my $default_qual = $$self{defaults}{QUAL};
    my $warn_sorted=1;
    my $warn_duplicates = exists($$self{warn_duplicates}) ? $$self{warn_duplicates} : 1;
    my $warn_total = 0;

    my ($prev_chrm,$prev_pos);

    for (my $i=0; $i<@$line; $i++)
    {
        if (!defined($$line[$i]) or $$line[$i] eq '' ) 
        {
            my $colname = $i<@{$$self{columns}} ? $$self{columns}[$i] : $i+1;
            $self->warn("The column $colname is empty at $$line[0]:$$line[1].\n");
            $warn->{EmptyColname} += 1;
            $warn->{line} += 1;
        }
    }

    my $x = $self->next_data_hash($line);
    if(validate_line($self, $x)){
        #moved this here to have all my warnings in the same spot
        $self->warn("Expected alphanumeric ID at $$x{CHROM}:$$x{POS}, but got [$$x{ID}]\n"); 
        $warn->{line} += 1;
        $warn->{MalformedChr} += 1;
    }

    # Is the position numeric?
    if ( !($$x{POS}=~/^\d+$/) ) { 
        $self->warn("Expected integer for the position at $$x{CHROM}:$$x{POS}\n");
        $warn->{line} += 1;
        $warn->{MalformedPos} += 1;
    }

    if ( $warn_duplicates ){
        if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos eq $$x{POS} ){
            $self->warn("Warning: Duplicate entries, for example $$x{CHROM}:$$x{POS}\n");
            $warn->{line} += 1;
            $warn->{Duplicate} += 1;
            $warn_duplicates = 0;
        }
    }

    # Is the file sorted?
    if ( $warn_sorted )
    {
        if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos > $$x{POS} ) 
        { 
            $self->warn("Warning: The file is not sorted, for example $$x{CHROM}:$$x{POS} comes after $prev_chrm:$prev_pos\n");
            $warn_sorted = 0;
            $warn->{line} += 1;
            $warn->{Sorting} += 1;
        }
        $prev_chrm = $$x{CHROM};
        $prev_pos  = $$x{POS};
    }

    # The reference base: one of A,C,G,T,N, non-empty.
    my $err = $self->validate_ref_field($$x{REF});
    if ( $err ){ 
        $self->warn("$$x{CHROM}:$$x{POS} .. $err\n");
        $warn->{line} += 1;
        $warn->{RefBase} += 1;
    }

    # The ALT field (alternate non-reference base)
    $err = $self->validate_alt_field($$x{ALT},$$x{REF});
    if ( $err ) {
        $self->warn("$$x{CHROM}:$$x{POS} .. $err\n");
        $warn->{line} += 1;
        $warn->{ALT_field} += 1;
    }

    # The QUAL field
    my $ret = $self->validate_float($$x{QUAL},$default_qual);
    if ( $ret ) { 
        $self->warn("QUAL field at $$x{CHROM}:$$x{POS} .. $ret\n");
        $warn->{line} += 1;
        $warn->{QUAL_field} += 1;
    }
    elsif ( $$x{QUAL}=~/^-?\d+$/ && $$x{QUAL}<-1 ) {
        $self->warn("QUAL field at $$x{CHROM}:$$x{POS} is negative .. $$x{QUAL}\n"); 
        $warn->{line} += 1;
        $warn->{QUAL_field} += 1;
    }

    # The FILTER field
    $err = $self->validate_filter_field($$x{FILTER});
    if ( $err ) {
        $self->warn("FILTER field at $$x{CHROM}:$$x{POS} .. $err\n");
        $warn->{line} += 1;
        $warn->{FILTER_field} += 1;
    }

    # The INFO field
    $err = $self->validate_info_field($$x{INFO},$$x{ALT});
    if ( $err ) {
        $self->warn("INFO field at $$x{CHROM}:$$x{POS} .. $err\n");
        $warn->{line} += 1;
        $warn->{INFO_field} += 1;
    } 

    #Genotype errors
    while (my ($gt,$data) = each %{$$x{gtypes}})
    {
        $err = $self->validate_gtype_field($data,$$x{ALT},$$x{FORMAT});
        if ( $err ) {
            $self->warn("column $gt at $$x{CHROM}:$$x{POS} .. $err\n");
            $warn->{line} +=1;
            $warn->{Genotypes} +=1;
        }
    }

    if ( scalar keys %{$$x{gtypes}} && (exists($$x{INFO}{AN}) || exists($$x{INFO}{AC})) )
    {
        my $nalt = scalar @{$$x{ALT}};
        if ( $nalt==1 && $$x{ALT}[0] eq '.' ) { $nalt=0; }
        my ($an,$ac) = $self->calc_an_ac($$x{gtypes},$nalt);    # Allow alleles in ALT which are absent in samples
        if ( exists($$x{INFO}{AN}) && $an ne $$x{INFO}{AN} ) 
        { 
            $self->warn("$$x{CHROM}:$$x{POS} .. AN is $$x{INFO}{AN}, should be $an\n"); 
            $warn->{line} +=1;
            $warn->{Genotypes} +=1;
        }
        if ( exists($$x{INFO}{AC}) && $ac ne $$x{INFO}{AC} ) 
        { 
            $self->warn("$$x{CHROM}:$$x{POS} .. AC is $$x{INFO}{AC}, should be $ac\n"); 
            $warn->{line} +=1;
            $warn->{Genotypes} +=1;
        }
    }

    return $warn;
}


sub validate_line{
    my ($self,$x) = @_;

    # Is the ID composed of alphanumeric chars
    if ( !($$x{ID}=~/^[\w;\.]+$/) ) { 
        return 1;
    }
    return 0;
}

