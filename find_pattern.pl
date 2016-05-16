#!/usr/bin/perl
# $Revision: 0.6 $
# $Date: 2015/03/14 $
# $Id: find_pattern.pl $
# $Author: Michael Bekaert $
# $Desc: Find fix allele patterns $
use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw/ uniq /;
use Data::Dumper;

#----------------------------------------------------------
our $VERSION = 0.6;
my $shift = 2;    #exported.haplotypes.tsv samples start at the third position!

#----------------------------------------------------------
sub reversed
{
    my $string = shift;
    if (defined $string)
    {
        $string = reverse(uc $string);
        $string =~ tr/ABCDGHMNRSTUVWXY/TVGHCDKNYSAABWXR/;
    }
    return $string;
}
my ($verbose, $fix, $refpop, $min, $ploidy, $min_snp, $max_snp, $grouping, $population, $haplofile, $whitefile, $arff, $genepop, $ade, $fasta, $mapfile, $snpfile, $vcfile) = (0, 0, 0, 0, 2, 1, 2, 0);
GetOptions(
           'haplotypes=s' => \$haplofile,
           'population=s' => \$population,
           'tag:s'        => \$mapfile,
           'snp:s'        => \$snpfile,
           'whitelist:s'  => \$whitefile,
           'min:i'        => \$min,
           'group:i'      => \$grouping,
           'arff:s'       => \$arff,
           'fasta:s'      => \$fasta,
           'genepop:s'    => \$genepop,
           'ade:s'        => \$ade,
           'minsnp:i'     => \$min_snp,
           'maxsnp:i'     => \$max_snp,
           'vcf:s'        => \$vcfile,
           'ref:i'        => \$refpop,
           'f|fix+'       => \$fix,
           'v|verbose!'   => \$verbose
          );
if (defined $haplofile && -r $haplofile && defined $population && -r $population && $min >= 0 && $grouping >= 0 && $grouping <= 2 && $refpop >= 0)
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
    my %vcf;
    if (defined $vcfile && -r $vcfile && open my $IN, q{<}, $vcfile)
    {
        while (<$IN>)
        {
            chomp;
            my @tmp = split m/\t/x;
            @{$vcf{$tmp[0]}} = ($tmp[2], $tmp[3], $tmp[1]) if (scalar @tmp >= 3 && substr($tmp[0], 0, 1) ne '#');
        }
        print {*STDERR} (scalar keys %vcf), " chromosome in the map!\n" if ($verbose);
        close $IN;
    }
    my (%physmap, %snpsmap);
    if (defined $mapfile && -r $mapfile && open my $IN, q{<}, $mapfile)
    {
        while (<$IN>)
        {
            chomp;
            my @tmp = split m/\t/x;
            if (scalar @tmp >= 10 && (!%whitelist || exists $whitelist{$tmp[2]}) && defined $tmp[9] && length $tmp[9] > 0)
            {
                if (defined $tmp[3] && length $tmp[3] > 0 && defined $tmp[4] && length $tmp[4] > 0 && defined $tmp[5] && length $tmp[5] > 0)
                {
                    $tmp[3] = $1 if ($tmp[3] =~ m/^(.*):\d+\.\.\d+$/g);
                    @{$physmap{$tmp[2]}} = ($tmp[9], $tmp[3], $tmp[4], $tmp[5]);
                }
                else { @{$physmap{$tmp[2]}} = $tmp[9]; }
            }
        }
        print {*STDERR} (scalar keys %physmap), " markers mapped!\n" if ($verbose);
        close $IN;
        if (defined $snpfile && -r $snpfile && open $IN, q{<}, $snpfile)
        {
            my ($nb, $last) = (0);
            while (<$IN>)
            {
                chomp;
                my @tmp = split m/\t/x;
                if (scalar @tmp >= 6 && (!%whitelist || exists $whitelist{$tmp[2]}))
                {
                    my @seq;
                    for my $i (5 .. (scalar @tmp - 1)) { push @seq, $tmp[$i] if (length($tmp[$i]) > 0); }
                    if (!defined $last || $last != $tmp[2]) { $nb = 0; }
                    else                                    { $nb++; }
                    @{$snpsmap{$tmp[2] . q{_} . chr(65 + $nb)}} = ($tmp[3], join(q{}, sort @seq));
                    $last = $tmp[2];
                }
            }
            print {*STDERR} (scalar keys %snpsmap), " SNP identified!\n" if ($verbose);
            close $IN;
        }
    }
    my (%pop, %group, %class);
    if (open my $IN, q{<}, $population)
    {
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
            foreach my $item (sort keys %class) { print {*STDERR} 'Group ', $item, ' [', $class{$item}, '] have ', $group{$class{$item}}, ' member', ($group{$class{$item}} > 1 ? q{s} : q{}), "\n"; }
        }
    }
    $refpop = 0 if (!exists $group{$refpop});
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
            #if (scalar @ind == scalar keys %pop)
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
                            last if ($item eq 'consensus' || $num_allele > $ploidy || $num_snp > $max_snp);
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
                                    foreach my $subitem (@{$alleles{$item}}) { push @{$all_alleles{$item}}, join(q{}, uniq(sort @{$subitem})); }
                                }
                                else
                                {
                                    for my $i (1 .. $num_snp) { push @{$all_alleles{$item}}, 'N'; }
                                }
                            }
                            undef %alleles;
                        }
                        if (%all_alleles && $num_snp >= $min_snp)
                        {
                            my $lasti = -1;
                            for my $i (0 .. ($num_snp - 1))
                            {
                                my $flag_fix = 1;
                                if ($fix)
                                {
                                    my @thelist;
                                    foreach my $item (keys %pop)
                                    {
                                        if (exists $all_alleles{$item}[$i])
                                        {
                                            push @thelist, $all_alleles{$item}[$i] if ($fix == 2 && $all_alleles{$item}[$i] ne 'N');
                                            $flag_fix = 0 if (length($all_alleles{$item}[$i]) > 1);
                                        }
                                    }
                                    $flag_fix = 0 if (@thelist && (scalar uniq sort @thelist) > 2);
                                }
                                if ($flag_fix)
                                {
                                    if ($grouping == 1 || $grouping == 2)
                                    {
                                        foreach my $trait (@traits)
                                        {
                                            my $true = 0;
                                            foreach my $refs (keys %pop)
                                            {
                                                $true = 0;
                                                my $ref = $all_alleles{$refs}[$i];
                                                next if ($ref eq 'N');
                                                $ref = $ref x $ploidy if (length($ref) == 1);
                                                if (length($ref) == $ploidy)
                                                {
                                                    my (%tmpline_arff, %tmpline_ade, %tmpline_genepop, %tmpline_fasta);
                                                    my $flag = $trait->{$refs};
                                                    foreach my $item (keys %pop)
                                                    {
                                                        if (exists $all_alleles{$item}[$i])
                                                        {
                                                            my $tmp = $all_alleles{$item}[$i];
                                                            $tmp = $tmp x $ploidy if (length($tmp) == 1);
                                                            if (length($tmp) == $ploidy)
                                                            {
                                                                if ($tmp eq 'NN' || ($flag eq $trait->{$item} && $ref eq $tmp) || ($flag ne $trait->{$item} && $ref ne $tmp)) { $true++; }
                                                                if (defined $arff)
                                                                {
                                                                    my $tmp2 = $tmp;
                                                                    $tmp2 =~ s/N+/\?/g;
                                                                    push @{$tmpline_arff{$item}}, $tmp2;
                                                                    push @{$header_arff{$data[0] . q{_} . chr(65 + $i)}}, $tmp2 if ($tmp2 ne '?');
                                                                }
                                                                if (defined $ade)
                                                                {
                                                                    my $tmp2 = $tmp;
                                                                    $tmp2 =~ tr/ATCGN/1234 /;
                                                                    push @{$tmpline_ade{$item}}, $tmp2;
                                                                }
                                                                if (defined $genepop)
                                                                {
                                                                    my $tmp2 = $tmp;
                                                                    $tmp2 =~ s/A/01/g;
                                                                    $tmp2 =~ s/C/02/g;
                                                                    $tmp2 =~ s/G/03/g;
                                                                    $tmp2 =~ s/T/04/g;
                                                                    $tmp2 =~ s/N/00/g;
                                                                    push @{$tmpline_genepop{$item}}, $tmp2;
                                                                }
                                                                if (defined $fasta) { push @{$tmpline_fasta{$item}}, $tmp; }
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
                                                        if (defined $arff) { delete $header_arff{$data[0] . q{_} . chr(65 + $i)}; }
                                                        if (!(defined $arff || defined $genepop || defined $fasta || defined $ade)) { delete $header_good{$data[0] . q{_} . chr(65 + $i)}; }
                                                    }
                                                    else
                                                    {
                                                        if ((defined $arff || defined $genepop || defined $fasta || defined $ade))
                                                        {
                                                            foreach my $item (keys %pop)
                                                            {
                                                                if (defined $arff) { push @{$line_arff{$item}}, @{$tmpline_arff{$item}}; }
                                                                if (defined $ade)  { push @{$line_ade{$item}},  @{$tmpline_ade{$item}}; }
                                                                if (defined $genepop) { push @{$line_genepop{$pop{$item}}{$item}}, @{$tmpline_genepop{$item}}; }
                                                                if (defined $fasta) { push @{$line_fasta{$item}}, @{$tmpline_fasta{$item}}; }
                                                            }
                                                        }
                                                        push @header,       $data[0] . q{_} . chr(65 + $i);
                                                        push @good_markers, $data[0];
                                                        last;
                                                    }
                                                }
                                            }
                                            last if ($true == (scalar keys %pop));
                                        }
                                    }
                                    else
                                    {
                                        my ($ref, $flag);
                                        foreach my $item (keys %pop)
                                        {
                                            if (exists $all_alleles{$item}[$i])
                                            {
                                                my $tmp = $all_alleles{$item}[$i];
                                                $tmp = $tmp x $ploidy if (length($tmp) == 1);
                                                if (length($tmp) == $ploidy)
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
                                                else
                                                {
                                                    print {*STDERR} $data[0], ': Ambibuity detected: ploidy unresolved!', "\n" if ($verbose);
                                                    undef @header;
                                                    last;
                                                }
                                            }
                                            $lasti = $i;
                                        }
                                    }
                                }
                            }
                        }
                    }
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
                    my $counter = 0;    ##
                    foreach my $item (keys %line_ade)
                    {
                        my $tmp = $item;
                        $tmp =~ s/[^a-zA-Z0-9]+//gx;

                        #print {*STDOUT} substr($tmp, 0, 8), "\t", join("\t", @{$line_ade{$item}}), "\n";
                        print {*STDOUT} substr($tmp, 0, 5), $counter++, "\t", join("\t", @{$line_ade{$item}}), "\n";    ##
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
                    if ($fix == 2 && scalar keys %vcf > 0)
                    {
                        print {*STDOUT} "##fileformat=VCFv4.1\n##handle=\n##batch=\n##reference=\n";
                        print {*STDOUT}
                          '##INFO=<ID=VRT,Number=1,Type=Integer,Description="Variation type, 1 - SNV: single nucleotide variation, 2 - DIV: deletion/insertion variation, 3 - HETEROZYGOUS: variable, but undefined at nucleotide level, 4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element, 6 - NO VARIATON: sequence scanned for variation, but none observed, 7 - MIXED: cluster contains submissions from 2 or more allelic classes, 8 - MNV: multiple nucleotide variation with alleles of common length greater than 1, 9 – Exception">',
                          "\n";
                        print {*STDOUT} '##INFO=<ID=FLANK-5,Number=1,Type=String,Description="5\' flanking sequence surrounding the variation.">', "\n";
                        print {*STDOUT} '##INFO=<ID=FLANK-3,Number=1,Type=String,Description="3\' flanking sequence surrounding the variation.">', "\n";
                        print {*STDOUT} '##INFO=<ID=CMT=1,Number=.,Type=String,Description="Comment.">',                                           "\n";
                        print {*STDOUT} "##FORMAT=<ID=NA,Number=1,Type=Integer,Description=\"Number of alleles for the population.\">\n";
                        print {*STDOUT} "##FORMAT=<ID=AC,Number=.,Type=Integer,Description=\"Allele count for each alternate allele.\">\n";
                        foreach my $pop (sort keys %class) { print {*STDOUT} "##population_id=", $pop, "\n"; }
                        print {*STDOUT} "#CHROM    	POS    	ID    	REF    	ALT    	QUAL	FILTER	INFO	FORMAT";
                        foreach my $pop (sort keys %class) { print {*STDOUT} "\t", $pop; }
                        print {*STDOUT} "\n";

                        foreach my $item (@header)
                        {
                            if (exists $header_good{$item})
                            {
                                my $id = substr $item, 0, -2;
                                if (%snpsmap && exists $snpsmap{$item} && %physmap && exists $physmap{$id} && exists $vcf{$physmap{$id}[1]})
                                {
                                    if (exists $header_good{$item}{$refpop})
                                    {
                                        my $ref_allele = substr(join(q{}, uniq(sort(@{$header_good{$item}{$refpop}}))), 0, 1);

                                        #my $ref_allele;
                                        #if (exists $header_good{$item}{$refpop}) {
                                        #                                    $ref_allele=substr(join(q{},uniq(sort(@{$header_good{$item}{$refpop}}))),0,1);
                                        #} else {
                                        # print {*STDERR} "Ignore '$item' as the reference allele is unknown!\n";
                                        #                                   $ref_allele=substr(join(q{},uniq(sort(@{$header_good{$item}{$refpop+1}}))),0,1);
                                        #}
                                        my $alt_allele = (substr($snpsmap{$item}[1], 0, 1) ne $ref_allele ? substr($snpsmap{$item}[1], 0, 1) : substr($snpsmap{$item}[1], 1, 1));
                                        if ($physmap{$id}[3] eq '+')
                                        {
                                            print {*STDOUT} $vcf{$physmap{$id}[1]}[2], "\t", ($physmap{$id}[2] + $snpsmap{$item}[0] + $vcf{$physmap{$id}[1]}[1]), "\t", $vcf{$physmap{$id}[1]}[0], ':g.',
                                              ($physmap{$id}[2] + $snpsmap{$item}[0] + $vcf{$physmap{$id}[1]}[1]), $ref_allele, q{>}, $alt_allele, "\t";
                                            print {*STDOUT} $ref_allele, "\t", $alt_allele, "\t.\tPASS\tVRT=1;FLANK-5=", substr($physmap{$id}[0], 0, $snpsmap{$item}[0]), ";FLANK-3=", substr($physmap{$id}[0], $snpsmap{$item}[0] + 1);
                                        }
                                        else
                                        {
                                            print {*STDOUT} $vcf{$physmap{$id}[1]}[2], "\t", ($physmap{$id}[2] - $snpsmap{$item}[0] + $vcf{$physmap{$id}[1]}[1]), "\t", $vcf{$physmap{$id}[1]}[0], ':g.',
                                              ($physmap{$id}[2] - $snpsmap{$item}[0] + $vcf{$physmap{$id}[1]}[1]), reversed($ref_allele), q{>}, reversed($alt_allele), "\t";
                                            print {*STDOUT} reversed($ref_allele), "\t", reversed($alt_allele), "\t.\tPASS\tVRT=1;FLANK-5=", reversed(substr($physmap{$id}[0], $snpsmap{$item}[0] + 1)), ";FLANK-3=",
                                              reversed(substr($physmap{$id}[0], 0, $snpsmap{$item}[0]));
                                        }
                                        print {*STDOUT} ';CMT="', $item, "\"\tNA:AC";
                                        foreach my $pop (sort keys %class)
                                        {
                                            print {*STDOUT} "\t";
                                            if (exists $header_good{$item}{$class{$pop}})
                                            {
                                                print {*STDOUT} (scalar @{$header_good{$item}{$class{$pop}}});
                                                print {*STDOUT} q{:};
                                                if (substr($header_good{$item}{$class{$pop}}[0], 0, 1) ne $ref_allele) { print {*STDOUT} (scalar @{$header_good{$item}{$class{$pop}}}); }
                                                else                                                                   { print {*STDOUT} '0'; }
                                            }
                                            else { print {*STDOUT} "0:0"; }
                                        }

                                        #                                    print {*STDOUT} "\tXXX" if (!exists $header_good{$item}{$refpop});
                                        print {*STDOUT} "\n";
                                    }
                                    else { print {*STDERR} "Ignore '$item' as the reference allele is unknown!\n"; }
                                }
                                else { print {*STDERR} "Skip '$item' as the chromosome is unknown!\n"; }
                            }
                        }
                    }
                    else
                    {
                        print {*STDOUT} 'Marker_SNP';
                        foreach my $pop (sort keys %group) { print {*STDOUT} "\tAllele_", $pop; }
                        print {*STDOUT} "\n";
                        foreach my $item (@header)
                        {
                            if (exists $header_good{$item})
                            {
                                my $id = substr $item, 0, -2;
                                print {*STDOUT} $item;
                                foreach my $pop (sort keys %group) { print {*STDOUT} "\t", (exists $header_good{$item}{$pop} ? q{\{} . join(q{,}, uniq(sort(@{$header_good{$item}{$pop}}))) . q{\}} : '{NN}'); }
                                if (%snpsmap && exists $snpsmap{$item} && %physmap && exists $physmap{$id})
                                {
                                    print {*STDOUT} "\t", substr($physmap{$id}[0], 0, $snpsmap{$item}[0]), q{[}, $snpsmap{$item}[1], q{]}, substr($physmap{$id}[0], $snpsmap{$item}[0] + 1);
                                    print {*STDOUT} "\t", $physmap{$id}[1], "\t", $physmap{$id}[2], "\t", $physmap{$id}[3] if (exists $physmap{$id}[1]);
                                }
                                print {*STDOUT} "\n";
                            }
                        }
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
      "Usage: $0 --haplotypes <batch_<num>.haplotypes.tsv> --population <popmap.txt>\nDescription: Test for diagnostic alleles or patterns between populations\n\n--haplotypes <file>\n    Raw haplotype file, automaticaly generated by Stacks, and called\n    batch_<num>.haplotypes.tsv.\n\n--population <file>\n    Population file used by Stacks.\n\n--tag <file>\n    batch_<num>.catalog.tags.tsv, required for physical mapping and marker sequences.\n\n--snp <file>\n    batch_<num>.catalog.snps.tsv, required for marker sequences.\n\n--whitelist <file>\n    Text file with the list of marker to only consider.\n\n--min <integer>\n    Minimum number of sample sharing alleles. [default 75%]\n\n--group <file>\n    Grouping [default 0].\n      0 between individuals [all];\n      1 between populations [species];\n      2 between groups of population [group].\n\n--arff <file>\n    Output as an ARFF file format.\n\n--fasta <file>\n    Output as a FASTA file format (SNP only).\n\n--genepop <file>\n    Output as an genepop file format.\n\n--minsnp <integer>\n    Minimum number of SNP. [default 1]\n\n--maxsnp <integer>\n    Maximum number of SNP. [default 2]\n\n--fix\n    Force fixed alleles only.\n\n--fix --fix\n    Force ONE fixed fallele only.\n\n--verbose\n    Becomes very chatty.\n\n";
}

#./find_pattern.pl --haplotypes batch_2.haplotypes.tsv --population farmed.txt -v --group 2 -d
#./find_pattern.pl --haplotypes tilapia.all.ref/batch_4.haplotypes.tsv --population population.named.txt -v --group 1 -maxsnp 2 -fix -fix --tag tilapia.all.ref/batch_4.catalog.tags.tsv --snp tilapia.all.ref/batch_4.catalog.snps.tsv --vcf=vcf_map.txt --ref=0
