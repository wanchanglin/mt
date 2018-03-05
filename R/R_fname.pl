#!/usr/bin/perl
# lwc-17-02-2009: extract functions in R scripts and giving index.
# Note: This is my first perl script!

use warnings;
use strict;

# @ARGV = qw/data_misc.R/;
my $no = 0;

while (<>) {
    #if(/^[^\s].*?<- function/){
    if(/^.*?<- function/){
        $no++;
        s/(^.*?)<-.*/## ($no). $1/g;
        print $_;
    }
}

## lwc-27-08-2015: what's for in the following code segment?
if (0){
    while (<>) {
        if(s/(^[^\s].*?)<- function.*/## $1/){
            print $_;
        }
    }
}
## lwc-17-02-2009: How to eval $no inside s///?
