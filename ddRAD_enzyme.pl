#!/usr/bin/perl
#
# Copyright 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# ddRAD_enzyme.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ddRAD_enzyme.pl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v3
# along with ddRAD_enzyme.pl. If not, see <http://www.gnu.org/licenses/>.
#
use strict;
use warnings;
use Getopt::Long;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;
use Bio::PrimarySeq;
use Bio::SeqIO;
my %enzyme = ();
$enzyme{'SbfI'} =   Bio::Restriction::Enzyme->new(-enzyme => 'SbfI', -seq => 'CCTGCA^GG');
$enzyme{'PstI'} =   Bio::Restriction::Enzyme->new(-enzyme => 'PstI', -seq => 'CTGCA^G');
$enzyme{'NlaIII'} = Bio::Restriction::Enzyme->new(-enzyme => 'NlaIII', -seq => 'CATG');
$enzyme{'SphI'} =   Bio::Restriction::Enzyme->new(-enzyme => 'SphI', -seq => 'GCATG^C');
$enzyme{'EcoRI'} =  Bio::Restriction::Enzyme->new(-enzyme => 'EcoRI', -seq => 'G^AATTC');
my ($verbose, $rad, $one, $two, $input) = (0, 0);
GetOptions(
           'i|input=s'  => \$input,
           '-1=s'       => \$one,
           '-2=s'       => \$two,
           'rad!'       => \$rad,
           'v|verbose!' => \$verbose
          );

if (   defined $input
    && -r $input
    && defined $one
    && exists $enzyme{$one}
    && defined $two
    && exists $enzyme{$two})
{
    my $enzymes = Bio::Restriction::EnzymeCollection->new(
                                    -enzymes => [$enzyme{$one}, $enzyme{$two}]);
    my $seq_obj = Bio::SeqIO->new(-file => $input, -format => 'fasta');
    my %total;
    my $skip = 0;
    my @all_frag;
    map { $total{$_->name} += 0 } $enzymes->each_enzyme;
    while (1)
    {
        my $seq = $seq_obj->next_seq || last;
        $skip++;
        my $ra =
          Bio::Restriction::Analysis->new(-seq => $seq, -enzymes => $enzymes);
        if ($rad)
        {
            print {*STDERR} $seq->id, "\n" if ($verbose);
            foreach my $i ($ra->fragment_maps($two))
            {
                my $sequence = Bio::PrimarySeq->new(-seq => $i->{seq});
                my $ra2 =
                  Bio::Restriction::Analysis->new(-seq     => $sequence,
                                                  -enzymes => $enzymes);
                if ($ra2->cuts_by_enzyme($one) > 0)
                {
                    my ($first, $last);
                    foreach my $j ($ra2->sizes($one))
                    {
                        $first = $j if (!defined $first);
                        $last = $j;
                    }
                    print {*STDOUT} $first, "\n", $last, "\n";
                }
            }
        }
        else
        {
            my $cuts = $ra->cut('multiple', $enzymes);
            print {*STDERR} $seq->id,
              '-> ', join(q{ }, $cuts->sizes('multiple_digest')), "\n"
              if ($verbose);
            map { print {*STDOUT} $_, "\n"; } $cuts->sizes('multiple_digest');
            map { push @all_frag, $_ } $cuts->sizes('multiple_digest')
              if ($verbose);
        }
        my $all_cutters = $ra->cutters;
        map { $total{$_->name} += $ra->cuts_by_enzyme($_->name) }
          $all_cutters->each_enzyme;
    }
    printf {*STDERR} "\n%-10s%s\n", 'Enzyme', 'Number of Cuts';
    foreach my $key (keys %total) {
        printf {*STDERR} "%-10s%d\n", $key, $total{$key};
    }
    print {*STDERR} (scalar @all_frag), " fragments\n\n" if ($verbose);
    my $name = ($rad ? 'rad' : 'cuts') . q{.} . $one . q{_} . $two;
    print {*STDERR}
      "#enzyme_cut(\"$name.txt\",450,550,\"$name.small.png\",\"$name.png\")\n\n";
}
else
{
    print {*STDOUT}
      "Usage $0 -i <input fasta genome> -1 <enzyme 1> -2 <enzyme2> [--rad] > list_of_fragment_size.txt\n\nOnly five restriction enzymes are encoded\n  SbfI:   CCTGCA^GG\n  PstI:   CTGCA^G\n  NlaIII: CATG^\n  SphI:   GCATG^C\n  EcoRI:  G^AATTC\n\n";
    print {*STDOUT}
      "#R commands:\nlibrary(ggplot2);\nenzyme_cut <- function(inputfile, sizemin=450, sizemax=550, smallpngfile=NULL,fullpngfile=NULL) {\n dist <-read.table(inputfile, col.names=c('val'), header=FALSE);\n dist <-transform(dist[grep(\"^\\\\d+\$\", dist\$val),,drop=F], val= as.numeric(as.character(val)));\n print(length(dist\$val));\n print(mean(dist\$val));\n print(median(dist\$val));\n dist2<-subset(dist, val <=sizemax);\n dist2<-subset(dist2, val >=sizemin);\n print(length(dist2\$val));\n if(!is.null(smallpngfile)) {\n   png(smallpngfile);\n   print(ggplot(dist, aes(x=val)) + stat_bin(binwidth=1) + scale_x_sqrt(limits=c(0,800)) + geom_vline(xintercept=sizemax)+geom_vline(xintercept=sizemin)+theme_bw());\n   dev.off();\n }\n if(!is.null(fullpngfile)) {\n   png(fullpngfile);\n   print(ggplot(dist, aes(x=val)) + stat_bin(binwidth=1) + scale_x_sqrt(limits=c(0,20000)) + geom_vline(aes(xintercept=median(val)),colour=\"#ec1c24\")+theme_bw());\n   dev.off();\n }\n}\nenzyme_cut(\"input.txt\",450,550,\"rad.small.png\",\"rad.png\")\n\n";
}
