#!/usr/bin/perl

use strict;
use LWP::Simple;

# Download Assembly records that are matching the input Enrez query
my $query = shift(@ARGV);
my $db = 'assembly';

#assemble the esearch URL
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

#post the esearch URL
my $output = get($url);
print "$url\n";

#parse WebEnv and QueryKey
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for ESearch-ESummary
#assemble the esummary URL
my $url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";

#post the esummary URL
my $docsums = get($url);
#~ print "$docsums";

### include this code for ESearch-EFetch
#assemble the efetch URL
my $url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web&retmode=text";
print "$url\n";

#post the efetch URL
#~ my $data = get($url);
#~ print "$data";
