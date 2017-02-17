#!/usr/bin/perl
#
# Copyright 2017, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# needle_matrix.pl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# needle_matrix.pl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v3
# along with needle_matrix.pl. If not, see <http://www.gnu.org/licenses/>.
#
use strict;
use warnings;
use Getopt::Long;
our $VERSION = 1.0;
my ($similarity, $keep, $gapopen, $gapextend, $fasta, $reference) = (0, 0, 10, 0.5);
GetOptions('i|in|fasta=s' => \$fasta, 'r|db|reference=s' => \$reference, 'k|keep!' => \$keep, 's|similarity!' => \$similarity, 'o|gapopen:f' => \$gapopen, 'e|gapextend:f' => \$gapextend);
if (defined $fasta && -r $fasta && defined $reference && -r $reference && defined $gapopen && defined $gapextend)
{
    use Bio::SeqIO;
    my $searched = ($similarity ? 'Similarity' : 'Identity');
    my $in = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');
    while (my $seq = $in->next_seq())
    {
        my $name = $seq->id();
        print {*STDERR} 'Processing: ', $name, '...';
        $name =~ s/\W+//gxm;
        my $out = Bio::SeqIO->new(-file => q{>} . $name . '.fa', -format => 'fasta');
        $out->write_seq($seq);
        $out->close();
        system('needle -asequence ' . $name . '.fa -bsequence ' . $reference . ' -outfile ' . $name . '.algn -brief -gapopen ' . $gapopen . ' -gapextend ' . $gapextend . ' >/dev/null 2>/dev/null') unless (-r $name . '.algn');
        if (-r $name . '.algn' && open(my $file, q{<}, $name . '.algn'))
        {
            print {*STDOUT} $name;
            while (<$file>)
            {
                if (/$searched.*\((\d+\.\d+)\%\)/oxm) { print {*STDOUT} "\t", $1, q{%}; }
            }
            print {*STDOUT} "\n";
            close $file;
            unlink $name . '.fa';
            unlink $name . '.algn' unless ($keep);
            print {*STDERR} " done\n";
        }
        else { print {*STDERR} " EMBOSS/Needle failed!\n"; }
    }
}
else
{
    print {*STDERR}
      "$0 release $VERSION\n\nRequiement: EMBOSS tools and BioPerl need to be installed\n e.g. apt-get install -y emboss bioperl --no-install-recommends\n\nUsage $0 -in <fasta sequences to align> -db <fasta sequences to align against> [...] > output.tsv\n\n--in <path>\n      Path th FASTA formated file with sequences to align. [mandatory]\n--db <path>\n      Path th FASTA formated file with sequences to align against. Can be the same as\n      the --in parametre is full matrix of pairwise comparison is needed. [mandatory]\n--gapopen <value>\n      NEEDLE --gapopen parametre. [default 10]\n--gapextend <value>\n      NEEDLE --gapextend parametre. [default 0.5]\n--keep\n      Keep the alignemnt file(s)\n--similarity\n      Report Similarity values rather than Identity (for protein only)\n\n";
}

#=======================================
#
# Aligned_sequences: 2
# 1: 1
# 2: 34
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 2321
# Identity:     899/2321 (38.7%)
# Similarity:  1182/2321 (50.9%)
# Gaps:         660/2321 (28.4%)
# Score: 3977.0
#
#
#=======================================
