#!/usr/bin/perl
# $Revision: 0.0 $
# $Date: 2014/02/14 $
# $Id: find_pattern.pl $
# $Author: Michael Bekaert $
# $Desc: Find fix allele patterns $
use strict;
use warnings;
use Getopt::Long;

#----------------------------------------------------------
our $VERSION = 0.0;
my $shift = 2;    #exported.haplotypes.tsv samples start at the third position!

#----------------------------------------------------------
my (
    $verbose,   $fix,      $min,  $threshold, $population, $haplofile,
    $whitefile, $fst_file, $arff, $genepop,   $fasta
   )
  = (0, 0, 0, 0.05);
GetOptions(
           'h|haplotypes=s' => \$haplofile,
           'p|population=s' => \$population,
           'w|whitelist:s'  => \$whitefile,
           'm|min:i'        => \$min,
           'fst:s'          => \$fst_file,
           't|threshold:f'  => \$threshold,
           'f|fix!'         => \$fix,
           'arff:s'         => \$arff,
           'fasta:s'        => \$fasta,
           'genepop:s'      => \$genepop,
           'v|verbose!'     => \$verbose
          );
if (   defined $haplofile
    && -r $haplofile
    && defined $population
    && -r $population
    && $min >= 0)
{
    my %whitelist;
    if (defined $whitefile && -r $whitefile && open my $IN, q{<}, $whitefile)
    {
        while (<$IN>)
        {
            chomp;
            my @tmp = split m/\t/x;
            $whitelist{$tmp[0]} = $tmp[0] if (exists $tmp[0]);
        }
        print {*STDERR} (scalar keys %whitelist), " markers in the whitelist!\n"
          if ($verbose);
        close $IN;
    }
    my %fst;
    if (defined $fst_file && -r $fst_file && (open my $IN, q{<}, $fst_file))
    {
        <$IN>;
        while (<$IN>)
        {
            chomp;
            my @tmp = split m/\t/x;

            # 6  'Column'
            # 9  'Fisher's P'
            # 17 'Corrected AMOVA Fst'
            if (   scalar @tmp >= 18
                && $tmp[9] <= $threshold
                && (!%whitelist || exists $whitelist{$tmp[1]}))
            {
                $fst{$tmp[1]}{$tmp[6]} = $tmp[17];
            }
        }
        close $IN;
        print {*STDERR} (scalar keys %fst), " Fsts collected!\n" if ($verbose);
    }
    my (%pop, %group, %weka);
    if (open my $IN, q{<}, $population)
    {
        my %class;
        while (<$IN>)
        {
            chomp;
            my @tmp = split m/\t/x;
            if (scalar @tmp >= 2 && $tmp[1] ne q{-})
            {
                if (!exists $class{$tmp[1]})
                {
                    $class{$tmp[1]} = scalar keys %class;
                    $group{$class{$tmp[1]}} = 0;
                }
                $pop{$tmp[0]}  = $class{$tmp[1]};
                $weka{$tmp[0]} = $tmp[0] . q{,} . $tmp[1];
                $group{$class{$tmp[1]}}++;
            }
        }
        close $IN;
        if ($verbose)
        {
            print {*STDERR} (scalar keys %pop), " samples to be used!\n";
            foreach my $item (keys %class)
            {
                print {*STDERR} 'Group ', $item, ' [', $class{$item},
                  '] have ', $group{$class{$item}}, " member\n";
            }
        }
    }
    $min = int((scalar keys %pop) * 0.75)
      if ($min == 0 || $min > scalar keys %pop);
    print {*STDERR} 'Threshold fixed at ', $min, "\n" if ($verbose);
}
else
{
    print
      "Usage: $0 --haplotypes <exported.haplotypes.tsv> --population <popmap.txt>\nDescription: Test for fixed allele or patterns between populations\n\n
--haplotypes <file>\n    Exported haplotypes by Stacks' export_sql.pl.\n\n
--population <file>\n    Population file used by Stacks.\n\n
--whitelist <file>\n    Text file with the list of marker to only consider.\n\n
--min <integer>\n    Minimum number of sample sharing alleles. [default 75%]\n\n
--fst <file>\n    \n\n
--threshold <float>\n    \n\n
--fix\n    Force fixed alleles only.\n\n
--arff <file>\n    \n\n
--fasta <file>\n    \n\n
--genepop <file>\n    \n\n
--verbose\n    Becomes very chatty.\n\n";
}
