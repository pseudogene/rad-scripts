#!/usr/bin/perl
# $Revision: 0.2 $
# $Date: 2014/02/25 $
# $Id: find_pattern.pl $
# $Author: Michael Bekaert $
# $Desc: Find fix allele patterns $
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw/ uniq /;
use Data::Dumper;

#----------------------------------------------------------
our $VERSION = 0.2;
my $shift = 2;    #exported.haplotypes.tsv samples start at the third position!

#----------------------------------------------------------
my ($debug, $verbose, $fix, $min, $ploidy, $max_snp, $grouping, $population, $haplofile, $whitefile, $arff, $genepop, $ade, $fasta) = (0, 0, 0, 0, 2, 5, 0);
GetOptions('h|haplotypes=s' => \$haplofile, 'p|population=s' => \$population, 'w|whitelist:s' => \$whitefile, 'm|min:i' => \$min, 'f|fix!' => \$fix, 'group:i' => \$grouping, 'arff:s' => \$arff, 'fasta:s' => \$fasta, 'genepop:s' => \$genepop, 'ade:s' => \$ade, 'v|verbose!' => \$verbose, 'd|debug!' => \$debug);
if (defined $haplofile && -r $haplofile && defined $population && -r $population && $min >= 0 && $grouping >= 0 && $grouping <= 2)
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
        print {*STDERR} (scalar keys %whitelist), " markers in the whitelist!\n" if ($verbose);
        close $IN;
    }
    my (%pop, %group);
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
                $pop{$tmp[0]} = $class{$tmp[1]};
                $group{$class{$tmp[1]}}++;
            }
        }
        close $IN;
        if ($verbose)
        {
            print {*STDERR} (scalar keys %pop), " samples to be used!\n";
            foreach my $item (sort keys %class) { print {*STDERR} 'Group ', $item, ' [', $class{$item}, '] have ', $group{$class{$item}}, " member\n"; }
        }
    }
    $min = int((scalar keys %pop) * 0.75) if ($min == 0 || $min > scalar keys %pop);
    print {*STDERR} 'Threshold fixed at ', $min, "\n" if ($verbose);
    if (%pop && %group)
    {
        my @traits;
        if ($grouping == 1) { push @traits, {%pop}; }
        elsif ($grouping == 2)
        {
            my @groups = keys %group;
            my $size   = scalar @groups;
            for (my $i = 0 ; $i < 2**$size ; $i++)
            {
                my $str = sprintf("%*.*b", $size, $size, $i);
                my %combination;
                for (my $j = 0 ; $j < $size ; $j++)
                {
                    #if (substr($str, $j, 1)) { $combination{$groups[$j]} = $groups[$j]; }
                    if (substr($str, 0, 1) ne '1' && substr($str, $j, 1)) { $combination{$groups[$j]} = $groups[$j]; }
                }
                if (scalar keys %combination > 0 && scalar keys %combination < $size)
                {
                    my %tmp;
                    foreach my $item (keys %pop) { $tmp{$item} = (exists $combination{$pop{$item}} ? 1 : 0); }
                    push @traits, {%tmp};
                }
            }
        }
        my @ind;
        if (open my $IN, q{<}, $haplofile)
        {
            my ($nb_all, $nb_selected, $nb_good, $nb_good_snp, $index) = (0, 0, 0, 0, 0);
            foreach my $item (split m/\t/x, <$IN>)
            {
                chomp $item;
                $index++;
                if ($index > $shift) { push @ind, $item; }
            }
            if (scalar @ind == scalar keys %pop)
            {
                my (@header, @good_markers);
                my (%line_arff, %line_genepop, %line_ade, %line_fasta, %header_arff, %header_good);
                while (<$IN>)
                {
                    $nb_all++;
                    chomp;
                    my @data = split m/\t/x;
                    if (scalar @data == (scalar @ind + $shift) && $data[1] >= $min && (!%whitelist || exists $whitelist{$data[0]}))
                    {
                        $nb_selected++;
                        my ($num_allele, $num_snp) = (0, 0);
                        my %alleles;
                        $index = 0;
                        foreach my $item (@data)
                        {
                            last if ($item eq 'consensus' || $num_allele > $ploidy || $num_snp > $max_snp);    ### consensus!!!!
                            $index++;
                            if ($index > $shift && $item ne '-')
                            {
                                my @item2 = sort(split(m/\//x, $item));
                                $num_allele = scalar @item2 if ($num_allele < scalar @item2);
                                $num_snp = length($item2[0]);
                                foreach my $allele (@item2)
                                {
                                    my $i = 0;
                                    foreach my $subitem (split m//x, $allele) { push @{$alleles{$ind[$index - ($shift + 1)]}[$i++]}, $subitem; }
                                }
                            }
                        }
                        my %all_alleles;
                        if ($index == scalar @data)
                        {
                            foreach my $item (keys %pop)
                            {
                                if (exists $alleles{$item})
                                {
                                    foreach my $subitem (@{$alleles{$item}}) { push @{$all_alleles{$item}}, join(q{}, sort @{$subitem}); }
                                }
                                else
                                {
                                    for my $i (1 .. $num_snp) { push @{$all_alleles{$item}}, 'N'; }
                                }
                            }
                            undef %alleles;
                        }
                        if (%all_alleles)
                        {
                            my $lasti = -1;
                            for my $i (0 .. ($num_snp - 1))
                            {
                                my $flag_fix = 1;
                                if ($fix)
                                {
                                    foreach my $item (keys %pop)
                                    {
                                        if (exists $all_alleles{$item}[$i]) { $flag_fix = 0 if (length($all_alleles{$item}[$i]) > 1); }
                                    }
                                }
                                if ($flag_fix)
                                {
                                    if ($grouping == 2)    #group    2
                                    {
## 2 ##
                                        foreach my $trait (@traits)
                                        {
                                            #if ($debug) {
                                            #    foreach my $item (keys $trait) { print {*STDERR} $trait->{$item}; }
                                            #    print {*STDERR} "\n";
                                            #}
                                            my ($ref, $flag);
                                            my ($true) = (0);
                                            foreach my $item (keys %pop)
                                            {
                                                if (exists $all_alleles{$item}[$i])
                                                {
                                                    my $tmp = $all_alleles{$item}[$i];
                                                    $tmp = $tmp x $ploidy if (length($tmp) == 1);
                                                    if (length($tmp) == $ploidy)
                                                    {
                                                        if (!defined $ref)
                                                        {
                                                            $ref  = $tmp;
                                                            $flag = $trait->{$item};
                                                        }
                                                        if (($flag eq $trait->{$item} && $ref eq $tmp) || ($flag ne $trait->{$item} && $ref ne $tmp)) { $true++; }

                                                        # test arff genepop fasta ade
                                                        if (!(defined $arff || defined $genepop || defined $fasta || defined $ade))
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ s/N+/\?/g;
                                                            push @{$header_good{$data[0] . q{_} . chr(65 + $i)}{$pop{$item}}}, $tmp2 if ($tmp2 ne '?');
                                                        }
                                                    }
                                                }
                                            }
                                            if ($true != (scalar keys %pop))
                                            {
                                                #remove header
                                                undef $header_good{$data[0] . q{_} . chr(65 + $i)};
                                            }
                                            else { push @header, $data[0] . q{_} . chr(65 + $i); }
                                        }
## /2 ##
                                    }
                                    else
                                    {
                                        my ($ref, $flag);
                                        my ($true) = (0);
                                        foreach my $item (keys %pop)
                                        {
                                            if (exists $all_alleles{$item}[$i])
                                            {
                                                my $tmp = $all_alleles{$item}[$i];
                                                $tmp = $tmp x $ploidy if (length($tmp) == 1);
                                                if (length($tmp) == $ploidy)
                                                {
                                                    if ($grouping == 1)    #species  1
                                                    {
## 1 ##
##no combinaison! what if order is wrong (pattern)
                                                        push @header, $data[0] . q{_} . chr(65 + $i) if ($i != $lasti);
                                                        push @good_markers, $data[0] if ($i != $lasti);
                                                        if (!defined $ref)
                                                        {
                                                            $ref  = $tmp;
                                                            $flag = $traits[0]{$item};
                                                        }
                                                        if (($flag eq $traits[0]{$item} && $ref eq $tmp) || ($flag ne $traits[0]{$item} && $ref ne $tmp)) { $true++; }

                                                        # test arff genepop fasta ade
                                                        if (!(defined $arff || defined $genepop || defined $fasta || defined $ade))
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ s/N+/\?/g;
                                                            push @{$header_good{$data[0] . q{_} . chr(65 + $i)}{$pop{$item}}}, $tmp2 if ($tmp2 ne '?');
                                                        }
## /1 ##
                                                    }
                                                    else    #all      0
                                                    {
                                                        push @header, $data[0] . q{_} . chr(65 + $i) if ($i != $lasti);
                                                        push @good_markers, $data[0] if ($i != $lasti);
                                                        if (defined $arff)
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ s/N+/\?/g;
                                                            push @{$line_arff{$item}}, $tmp2;
                                                            push @{$header_arff{$data[0] . q{_} . chr(65 + $i)}}, $tmp2 if ($tmp2 ne '?');
                                                        }
                                                        if (defined $ade)
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ tr/ATCGN/1234 /;
                                                            push @{$line_ade{$item}}, $tmp2;
                                                        }
                                                        if (defined $genepop)
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ s/A/01/g;
                                                            $tmp2 =~ s/C/02/g;
                                                            $tmp2 =~ s/G/03/g;
                                                            $tmp2 =~ s/T/04/g;
                                                            $tmp2 =~ s/N/00/g;
                                                            push @{$line_genepop{$pop{$item}}{$item}}, $tmp2;
                                                        }
                                                        if (defined $fasta) { push @{$line_fasta{$item}}, $tmp; }
                                                        if (!(defined $arff || defined $genepop || defined $fasta || defined $ade))
                                                        {
                                                            my $tmp2 = $tmp;
                                                            $tmp2 =~ s/N+/\?/g;
                                                            push @{$header_good{$data[0] . q{_} . chr(65 + $i)}{$pop{$item}}}, $tmp2 if ($tmp2 ne '?');
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    print {*STDERR} $data[0], ': Ambibuity detected: ploidy unresolved!', "\n" if ($verbose);
                                                    undef @header;
                                                    last;
                                                }
                                            }
                                            $lasti = $i;
                                        }
                                        if ($grouping == 1 && $true != (scalar keys %pop))    #species  1
                                        {
## 1 ##
                                            #remove header
                                            undef $header_good{$data[0] . q{_} . chr(65 + $i)};
## /1 ##
                                        }
                                    }
                                }
                            }
                        }
                    }
                    last if ($nb_selected > 3 && $debug);                                     ##
                }
                $nb_good = scalar uniq sort @good_markers;
                undef @good_markers;
                $nb_good_snp = scalar @header;
                if (defined $arff && @header && scalar @header > 0)
                {
                    print {*STDOUT} "\@RELATION $haplofile\n\n\@ATTRIBUTE Sample STRING\n\@ATTRIBUTE Population {", join(q{,}, uniq(sort(keys %group))), "}\n";
                    foreach my $item (@header) { print {*STDOUT} '@ATTRIBUTE ', $item, ' {', join(q{,}, uniq(sort(@{$header_arff{$item}}))), "}\n"; }
                    print {*STDOUT} "\n\@DATA\n";
                    foreach my $item (keys %line_arff) { print {*STDOUT} $item, q{,}, $pop{$item}, q{,}, join(q{,}, @{$line_arff{$item}}), "\n"; }
                }
                if (defined $ade && @header && scalar @header > 0)
                {
                    print {*STDOUT} "Samples\t", join("\t", @header), "\n";
                    foreach my $item (keys %line_ade)
                    {
                        my $tmp = $item;
                        $tmp =~ s/[^a-zA-Z0-9]+//gx;
                        print {*STDOUT} substr($tmp, 0, 8), "\t", join("\t", @{$line_ade{$item}}), "\n";
                    }
                    my @tmp;
                    foreach my $item (keys %line_ade) { push @tmp, $pop{$item}; }
                    print {*STDERR} "R/adegenet population vector:\n pop <- c('", join('\',\'', @tmp), "');\n";
                }
                if (defined $genepop && @header && scalar @header > 0)
                {
                    print {*STDOUT} "Find_pattern.pl version $VERSION; Genepop 4.+\n";
                    print {*STDOUT} join(q{,}, @header), "\n";
                    foreach my $pop (keys %line_genepop)
                    {
                        print {*STDOUT} "pop\n";
                        foreach my $item (keys %{$line_genepop{$pop}}) { print {*STDOUT} $item, ",\t", join("\t", @{$line_genepop{$pop}{$item}}), "\n"; }
                    }
                }
                if (defined $fasta && @header && scalar @header > 0)
                {
                    foreach my $item (keys %line_fasta) { print {*STDOUT} q{>}, $item, ' [', $pop{$item}, "]\n", join(q{}, @{$line_fasta{$item}}), "\n"; }
                }
                if (!(defined $arff || defined $genepop || defined $fasta || defined $ade) && @header && scalar @header > 0)
                {
                    print {*STDOUT} 'Marker_SNP';
                    foreach my $pop (sort keys %group) { print {*STDOUT} ',Allele_', $pop; }
                    print {*STDOUT} "\n";
                    foreach my $item (@header)
                    {
                        print {*STDOUT} $item;
                        foreach my $pop (sort keys %{$header_good{$item}}) { print {*STDOUT} q{,}, (exists $header_good{$item}{$pop} ? q{\{} . join(q{,}, uniq(sort(@{$header_good{$item}{$pop}}))) . q{\}} : q{}); }
                        print {*STDOUT} "\n";
                    }
                }
            }
            close $IN;
            print {*STDERR} "Total markers read: $nb_all\nMarker analysed:    $nb_selected\nMarker selected:    $nb_good\nSNP selected:       $nb_good_snp\n\n";
        }
    }
}
else
{
    print
      "Usage: $0 --haplotypes <batch_<num>.haplotypes.tsv> --population <popmap.txt>\nDescription: Test for fixed alleles or patterns between populations\n\n--haplotypes <file>\n    Raw haplotype file, automaticaly generated by Stacks, and called\n    batch_<num>.haplotypes.tsv.\n\n--population <file>\n    Population file used by Stacks.\n\n--whitelist <file>\n    Text file with the list of marker to only consider.\n\n--min <integer>\n    Minimum number of sample sharing alleles. [default 75%]\n\n--group <file>\n    Grouping [default 2].\n      0 between individuals;\n      1 between populations;\n      2 between groups of population.\n\n--arff <file>\n    Output as an ARFF file format.\n\n--fasta <file>\n    Output as a FASTA file format (SNP only).\n\n--genepop <file>\n    Output as an genepop file format.\n\n--fix\n    Force fixed alleles only.\n\n--verbose\n    Becomes very chatty.\n\n";
}

#grouping
#old        new
#all        0        done!
#group      2
#species    1        in progress!
#./find_pattern.pl --haplotypes batch_2.haplotypes.tsv --population farmed.txt -v --group 1 -d
#./find_pattern.pl --haplotypes batch_2.haplotypes.tsv --population farmed.txt -v --group 2
