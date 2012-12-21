# Pneumotyper
Copyright Frédéric Raymond 2012

Pneumotyper is a software tool to interpreet results obtained from a Streptococcus pneumoniae serotyping assay that is currently being submitted to a peer-reviewed journal. This software takes as an input the results from the AutoGenomics INFINITI "SPNEUMO" assay.

The software requires Perl. 

## Command-line use 

perl Pneumotyper.pl Results.txt

## Results interpretation

Example : 
MA87708 2012092019072317        Pos(2.52)       Pos(48.28)      Pos(40.51)      0.5     19a

Column 1 : Sample ID (MA87708)
Column 2 : Sample Run ID (2012092019072317)
Column 3 : Internal control. Pos = Positive. Neg = Negative. Number in parenthesis indicates the observed signal.
Column 4 : Pneumolysin probe. Pos = Positive. Neg = Negative. Number in parenthesis indicates the observed signal.
Column 5 : Autolysin probe . Pos = Positive. Neg = Negative. Number in parenthesis indicates the observed signal.
Column 6 : Serotyping score
Column 7 : Most probable serotype.

## Licence

The Pneumotyper software distributed under the GPLv3 Licence.

 This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

                  
