#!/usr/bin/perl -w
use strict;
use warnings;

##    Copyright Frédéric Raymond 2012
##    The Pneumotyper software is distributed under the GPLv3 Licence. See the LICENCE file for full licence details.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>. 

##### Output variables.
## The value of $toprint allows to change the format of the output. Other formats than best are used for debugging purposes. 

my $toprint = "best";      ### all ou best
#my $toprint = "all";      ### all ou best
#my $toprint = "stats";
#my $toprint = "new";
#### Fonction pour lire les fichiers de resultats d'AutoGenomics

### Exemple of input file
#ID:2011012611264193,Sample:-GECHIP-R169-1
#Assay: SPNEUMO8
#Analyte                   Mean            Ratio          Analysis
#-----------------------------------------------------------------------
#SpneumoD003-B             18.62            0.12          neg
#SpneumoD004-B             14.61            0.10          neg
#...
#SpneumoD026-Pn          1539.69           10.25          POSITIVE
#SpneumoD046-Ly           633.64            4.22          POSITIVE
#IC1                       10.55            0.07          neg
#IC2                       14.41            0.10          neg
#HYBC                     381.09            0.00          
#BKGD                     150.14            0.00          
#----------------Test-Errors--------------------------------------

### Sorting subroutine

sub par_num ($$) {
        my ($gauche, $droit) = @_;
        return $gauche <=> $droit;
}

## The software requires the SerotypeDatabase.tab file to compute serotypes.

my %results;
my @id;
my @sample;
my $begin = "false";
my @positions = ("B", "C", "D", "E", "F", "I", "J", "M", "N", "P", "Q", "R");

my %means;

my $resultfile = shift;

#my %sort_score;
#my %sero_results;

my $sumlow;
my $countlow;

open (RESULTS, $resultfile) || die ("Could not open AG result file");
while(<RESULTS>){
        chomp $_;    
        if($_ =~ /^ID/){
                my @header = split(/\,/, $_);
                @id = split(/\:/, $header[0]);
                @sample = split(/\:/, $header[1]);
                $sample[1] =~ s/\t//g;
                $sample[1] =~ s/ //g;
                $sample[1] =~ s/^-//;
                $results{$id[1]}{"Sample"} = $sample[1];
                $results{$id[1]}{"None"} = 0.01;
                foreach my $position (@positions){
                         $results{$id[1]}{$position}="None";
                }
                $begin = "false";
                $sumlow=0;
                $countlow=0;
        }
        if($_ =~ /BKGD/){
                $begin = "false";
                if ($countlow > 0){
                       $means{$id[1]} = $sumlow/$countlow;
                } else {
                        $means{$id[1]} = 0.25;
                }
                if($toprint eq "stats"){
                        print "$sumlow\t$countlow\t$means{$id[1]}\n";
                }
        }
        if ($begin eq "true" && $_ !~ /^---------/){
                my @line = split(/\t/, $_);
                $line[0] =~ s/ //g;
                $line[2] =~ s/\t//g;
                $line[2] =~ s/ //g;
                $line[2] =~ s/,/./g;
                if($line[2]>0.01 && $line[2]<0.7){
                        $countlow++;
                        $sumlow = $sumlow + $line[2];
                }
        }
        if($_ =~ /Analyte/){
                if($begin eq "false"){
                        $begin = "true";
                }       
        }

}


open (RESULTS, $resultfile) || die ("Could not open AG result file");
while(<RESULTS>){
        chomp $_;
        if($_ =~ /^ID/){
                my @header = split(/\,/, $_);
                @id = split(/\:/, $header[0]);
                @sample = split(/\:/, $header[1]);
                $sample[1] =~ s/\t//g;
                $sample[1] =~ s/ //g;
                $sample[1] =~ s/^-//;
                $results{$id[1]}{"Sample"} = $sample[1];
                $results{$id[1]}{"None"} = 0;
                foreach my $position (@positions){
                         $results{$id[1]}{$position}="None";
                }
                $begin = "false";
        }
        my $buffer = $means{$id[1]};
        if($_ =~ /^Assay/){
                my @assay = split(/\:/, $_);
                $results{$id[1]}{"Assay"} = $assay[1];
        }
        if($_ =~ /BKGD/){
                $begin = "false";
        }
        if ($begin eq "true" && $_ !~ /^---------/){
                my @line = split(/\t/, $_);
                $line[0] =~ s/ //g;
                $line[2] =~ s/\t//g;
                $line[2] =~ s/ //g;
                $line[2] =~ s/,/./g;
                my @setgroup = split(/-/, $line[0]);
                $results{$id[1]}{$line[0]} = $line[2];
                if(exists($setgroup[1])){
                        if(exists($results{$id[1]}{$setgroup[1]})){
                                my $temp = $results{$id[1]}{$setgroup[1]};
                                my $oldscore = $results{$id[1]}{$temp};
                                if($line[2] > $oldscore){
                                        if($line[2] > $oldscore+$buffer){ 
                                                $results{$id[1]}{$setgroup[1]} = $line[0];
                                        } else {
                                                $results{$id[1]}{$setgroup[1]} = "$setgroup[1]Ambiguous";
                                                $results{$id[1]}{"$setgroup[1]Ambiguous"} = $line[2];
                                        }
                                }
                        }
                }
        }
        if($_ =~ /Analyte/){
                if($begin eq "false"){
                        $begin = "true";
                }
        }
}
close (RESULTS);

##### Load serotype definition

my $serotypefile = "SerotypeDatabase.tab";
my %serotypes;

open (SEROTYPES, $serotypefile) || die ("Could not open AG result file");
while(<SEROTYPES>){
        chomp $_;
        $_ =~ s/SpneumoD009-025-F/SpneumoD009025-F/;
        my @description = split(/\t/, $_);
        my $i = 0;
        my $count=0;
        foreach my $position (@positions){
                $i++;
                $serotypes{$description[0]}{$position}= $description[$i];
                if($description[$i] =~ /Spneumo/){
                        $count++;
                }
        }
        $serotypes{$description[0]}{Count} = $count;
}
close(SEROTYPES);

## Serotype determination

my @echantillons = keys %results;
my @sero = keys %serotypes;

my $buffer = $means{$id[1]};

my $threshold = 0.5;
my $threshold_low = $buffer;
my $thresholdpneumo = 1.5;

foreach my $echantillon (@echantillons){
        my $scoringmargin = "";
        my %sort_scores;  
        my %sero_results;
        my $toprintnewscore = 0;
        my %analysis;
        my $bestscore = 0;
        my $bestscore2 = 0;
        my $bestlist = "Negative";
        my $bestlist2 = "Negative";
        my $pne = "neg";
        my $lytA = "neg";
        my $IC = "neg";
        my $infos = "0\t0\t0\t0\t0\t0\t0\t0\t0";

        ### Corriger probe 009 025
     if($results{$echantillon}{F} eq "SpneumoD009-F" || $results{$echantillon}{F} eq "Spneumo025-F"){
                my $temp = $results{$echantillon}{F};
                $results{$echantillon}{F} = "SpneumoD009025-F";
                $results{$echantillon}{$results{$echantillon}{F}} = $results{$echantillon}{$temp};
        }

        ### Imprimer details de l'echantillon

        ### Verifier Pneumolysin
        if($results{$echantillon}{'SpneumoD026-Pn'}>=$thresholdpneumo){
                $pne = "Pos($results{$echantillon}{'SpneumoD026-Pn'})";
        } else {
                $pne = "neg($results{$echantillon}{'SpneumoD026-Pn'})";
        }

        ### Verifier Autolysin
        if($results{$echantillon}{'SpneumoD046-Ly'}>=$thresholdpneumo){
                $lytA = "Pos($results{$echantillon}{'SpneumoD046-Ly'})";
        } else {
                $lytA = "neg($results{$echantillon}{'SpneumoD046-Ly'})";
        }

        ### Verifier Controle Interne
        if($results{$echantillon}{'IC1'}>=$threshold){
                $IC = "Pos($results{$echantillon}{'IC1'})";
        } else {
                $IC = "neg($results{$echantillon}{'IC1'})";
        }

        ### Comparer avec les definintions de serotypes ;
        if($pne =~ /Pos/ || $lytA =~ /Pos/){
        foreach my $sero (@sero){
                my $concordant = 0;
                my $PN = 0;
                my $NP = 0;
                my $discordant = 0;
                my $concordant_low = 0;
                my $discordant_low = 0;
                my $newscore = 0;
                my $score = 0;
                my $NN = 0;
                foreach my $pos(@positions){
                        my $samplepos;
                        my $seropos;

                        if($results{$echantillon}{$results{$echantillon}{$pos}}>=$threshold){
                                $samplepos = 1;
                        } else {
                                $samplepos = 0;
                        }
                        if($serotypes{$sero}{$pos} =~ /Spneumo/){
                                $seropos = 1;
                        } else {
                                $seropos = 0;
                        }
                        if($samplepos == 1 && $seropos == 0){
                                $PN++;
                                if($results{$echantillon}{$results{$echantillon}{$pos}}<1){
                                        $newscore = $newscore - $results{$echantillon}{$results{$echantillon}{$pos}};
                                } else {            
                                        $newscore = $newscore - 1;
                                }
                        }
                        if($samplepos == 0 && $seropos == 1){
                                $NP++;
                        }
                       if($samplepos == 0 && $seropos == 0){
                                $NN++;
                        }
                        if($samplepos == 1 && $seropos == 1 && $results{$echantillon}{$pos} ne $serotypes{$sero}{$pos}){
                                $discordant++;
                                if($results{$echantillon}{$results{$echantillon}{$pos}}<1){
                                        $newscore = $newscore - $results{$echantillon}{$results{$echantillon}{$pos}};
                                } else {
                                        $newscore = $newscore - 1;
                                }
                        }
                        if($samplepos == 1 && $seropos == 1 && $results{$echantillon}{$pos} eq $serotypes{$sero}{$pos}){
                                if($results{$echantillon}{$results{$echantillon}{$pos}}>=0.1){
                                        $concordant++;
                                        if($results{$echantillon}{$results{$echantillon}{$pos}}<1){
                                                $newscore = $newscore + $results{$echantillon}{$results{$echantillon}{$pos}};
                                        } else {
                                                 $newscore = $newscore + 1;
                                        }
                                }
                        }
                        if($results{$echantillon}{$results{$echantillon}{$pos}}>=$threshold_low && $results{$echantillon}{$results{$echantillon}{$pos}} < $threshold && $seropos == 1 && $results{$echantillon}{$pos} eq
$serotypes{$sero}{$pos}){
                                $concordant_low++;
                                #$concordant_low = $concordant_low + $results{$echantillon}{$results{$echantillon}{$pos}};
                        }
                        if($results{$echantillon}{$results{$echantillon}{$pos}}>=$threshold_low && $results{$echantillon}{$results{$echantillon}{$pos}} < $threshold && $seropos == 1 && $results{$echantillon}{$pos} ne 
$serotypes{$sero}{$pos}){
                                $discordant_low++;
                        }
                        
                }

### Scoring

                  if(($serotypes{$sero}{Count}+$NP+$NN+$discordant+$discordant_low)>0){
                           $newscore = ($concordant + $concordant_low - $PN)/($serotypes{$sero}{Count}+$NP+$NN+$discordant+$discordant_low);
                 } else {
                        $newscore = 0;
                }
                if ($newscore<0){
                        $newscore = 0;
                }
                if($serotypes{$sero}{Count}-$NP>0){
                        $score=$newscore;
                        $toprintnewscore = ($concordant+$concordant_low - $PN)/($serotypes{$sero}{Count} + $discordant + $discordant_low+$NP);
                }

                ### To print all serotype results.
                if($toprint eq "all"){
                        if($pne =~ /Pos/){
                          print $results{$echantillon}{Sample};
                          print "\t$echantillon\t$sero\t$concordant\t$serotypes{$sero}{Count}\t$PN\t$NP\t$NN\t$discordant\t$concordant_low\t$discordant_low\t";
                          print "$score\t$toprintnewscore\t$buffer";
                          print "\n";
                       }
                }
                if(exists($sort_scores{$score})){
                        $sort_scores{$score} = $sort_scores{$score} . ";$sero";
                } else {
                        $sort_scores{$score} = "$sero";
                }
                $sero_results{$sero} = "$concordant\t$serotypes{$sero}{Count}\t$PN\t$NP\t$discordant\t$concordant_low\t$discordant_low\t$score\t$toprintnewscore\t$buffer";
          
                ### To print the best serotype result.
                if($toprint eq "best"){ 
                        my @serosplit = split(/_/, $sero);
                        if($score==$bestscore){
                                $bestlist = $bestlist . ";" . $serosplit[0];
                        }
                        if($score>$bestscore){
                                if($score>=0.4){
                                        $bestscore = $score;
                                        $bestlist = $serosplit[0];
                                        $infos = "$concordant\t$serotypes{$sero}{Count}\t$PN\t$NP\t$discordant\t$concordant_low\t$discordant_low\t$score\t$toprintnewscore\t$buffer";
                                }
                                if($score<0.4){
                                        $bestscore = $score;
                                        $bestlist = "Negative\/" . $serosplit[0];
                                        $infos = "$concordant\t$serotypes{$sero}{Count}\t$PN\t$NP\t$discordant\t$concordant_low\t$discordant_low\t$score\t$toprintnewscore\t$buffer";
                                }
                        }
                        if($toprintnewscore == $bestscore2){
                                $bestlist2 = $bestlist2 . ";" . $serosplit[0];
                        }
                        if($toprintnewscore > $bestscore2){
                                if($score>=0.4){
                                        $bestscore2 = $toprintnewscore;
                                        $bestlist2 = $serosplit[0];
                                }
                                if($score<0.4){
                                        $bestscore2 = $toprintnewscore;
                                        $bestlist2 = "Negative\/" . $serosplit[0];
                                }
                        }
                }
        } 
        if($bestscore>0){
                my @serokeyz = keys %sort_scores;
                my @sorted_sero = sort par_num @serokeyz;
                my $countnb = @serokeyz;
                my $difference =  $sorted_sero[$countnb-1] - $sorted_sero[$countnb-2];
                if($difference<0.01 && $countnb>0){
                        my @serosplit = split(/;/, $sort_scores{$sorted_sero[$countnb-2]});
                        foreach my $seroX (@serosplit){
                                if($seroX =~ /_/){
                                        my @splitit = split(/_/, $seroX);
                                        $scoringmargin = $scoringmargin . ";$splitit[0]";
                                } else {
                                        $scoringmargin = $scoringmargin . ";$seroX";
                                }
                        }
                }
        }
        }

        if(($pne =~ /Pos/ || $lytA =~ /Pos/) && $bestscore <= 0){
                $bestlist = "Untyped";
        }
        if(($pne =~ /Pos/ || $lytA =~ /Pos/) && $bestscore2 <= 0){
                $bestlist2 = "Untyped";
        }

### To print the best serotype result.
        if($toprint eq "best"){
                print "$results{$echantillon}{Sample}\t$echantillon\t";
                print "$IC\t$pne\t$lytA\t"; 
                if($bestlist eq $bestlist2){
                        print "$bestscore\t$bestlist";
                } else {
                        print "$bestscore\t$bestlist;$bestlist2";
                }
                if ($scoringmargin ne ""){
                        print "$scoringmargin";
                }
                print "\n"
        }
}
