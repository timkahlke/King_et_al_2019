#!/usr/bin/perl

use strict;
use Getopt::Std;
use Data::Dumper;


eval{
    _main();
};


if($@){
    print $@;
}


sub _main{
    our %opts;
    getopts("b:f:p:l:o:m:",\%opts);
    my $bf = $opts{'b'};
    my $ff = $opts{'f'};
    my $pi = $opts{'p'}?$opts{'p'}:0;
    my $ml = $opts{'l'}?$opts{'l'}:0;
    my $overlap = $opts{'m'}?$opts{'m'}:0;
    my $of = $opts{'o'};

    if(!$bf||!$ff||!$of){
        _usage();
    }

    my $hits = _read_blast($bf,$pi,$ml);
    open(my $oh,">",$of) or die "Failed to open $of";
    open(my $ih,"<",$ff) or die "Failed to open file $ff";


    my ($seqs,$name,$ind);
    while(my $line=<$ih>){
        if($line=~/^>.*$/){
            if($hits->{$name}){
                my $hit_seqs = _extract_hit_seq($hits->{$name},$seqs,$overlap);
                foreach my $s(keys(%$hit_seqs)){
                    print $oh ">$name.$s\n".$hit_seqs->{$s}."\n\n";
                }
            }
            my @split = grep{$_ ne ""}split(/[\s\t\n|]/,$line);

            # guess structure of database name string, e.g. emb|XFBCS.1|blablabla
            # will set ind = 1
            if(!(defined($ind))){
                for(my $i=0;$i<scalar(@split);$i++){
                    if($hits->{@split[$i]}){
                        $ind = $i;
                    }
                }
            }

            next unless defined($ind);
            $name = @split[$ind];
            $seqs = "";
        }
        else{
            $line=~s/[\s\t\n\r]//g;
            $seqs.=$line;
        }
    }

    if($hits->{$seqs}){
        my $hit_seqs = _extract_hit_seqs($hits->{$name},$seqs,$overlap);
        foreach my $s(keys(%$hit_seqs)){
            print $oh ">$name.$s\n".$hit_seqs->{$s}."\n\n";
        }
    }


    close($oh);
}


sub _extract_hit_seq{
    my ($hits,$seqs,$overlap) = @_;
    my $hit_seqs = {};

    foreach my $s(@$hits){
        my ($start,$stop);
        if($overlap){
            $start = _max($s->{'start'}-$overlap,1);
            $stop = _min($s->{'stop'}+$overlap,length($seqs->{$s}));
        }
        else{
            $start = $s->{'start'};
            $stop = $s->{'stop'};
        }
        $hit_seqs->{$start."_".$stop} = substr($seqs,$start-1,$stop-$start);
    }
    return $hit_seqs;
}


sub _read_blast{
    my ($file,$pi,$ml) = @_;
    my $hits = {};
    open(my $bh,"<",$file) or die "Failed to open $file";
    while(my $line=<$bh>){
        $line=~s/[\r\n]//g;
        my @s = split(/\t/,$line);
        if($pi){
            if($s[2]<$pi){;
                next;
            }
        }
        if($ml){
            if($s[3]<$ml){
                next;
            }
        }
        if($hits->{$s[1]}){
            push @{$hits->{$s[1]}},{'start'=>_min($s[8],$s[9]),'stop'=>_max($s[8],$s[9])};
        }
        else{
            $hits->{$s[1]} = [{'start'=>_min($s[8],$s[9]),'stop'=>_max($s[8],$s[9])}];
        }
    }
    close($bh);
    return $hits;
}

sub _max{
    my ($a,$b) = @_;
    if($a<$b){
        return $b;
    }
    return $a;
}

sub _min{
    my ($a,$b) = @_;
    if($a>$b){
        return $b;
    }
    return $a;
}

sub _usage{
    print STDOUT "\n\nScript to extract a fasta file from hits by returning just the hit region of each hit.\n";
    print STDOUT "Parameter:\n";
    print STDOUT "b : blast file\n";
    print STDOUT "l : (optional) minimum length of hit\n";
    print STDOUT "p : (optional) minimum percent id of hit\n";
    print STDOUT "m : (optional) overlap, i.e., additional upstream and downstream to be extracted\n";
    print STDOUT "f : fasta file of blast database\n";
    print STDOUT "o : output fasta file\n";
    print STDOUT "\n\n";
    exit;
}

