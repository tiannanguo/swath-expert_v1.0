#!/usr/bin/perl
#tiannan,2015,IMSB,ETH
use Cwd;

$in=shift; #"feature_alignment_requant_matrix.tsv";
$in2=shift; #"140306PC-DDA-92files_step4.tsv";
$o=shift ; #"140306PC-DDA-92files_step5.tsv";

my %tg;
@d=&oF($in);
$title=shift @d;
foreach my $d(@d){
          my @dd=&s($d);
          if ($dd[0]=~/^(\d+_.*)\_run\d/){
               $tg{$1}=1;
          }
}

@d2=&oF($in2);
$title2=shift @d2;
open(OUT,">$o");
print OUT "$title2\n";
foreach my $d(@d2){
          my @dd=&s($d);
          if ($tg{$dd[6]}){
               print OUT "$d\n";
          }
}
close OUT;

sub oF{
    my $file=shift;
    my @d=read_file($file);
    foreach (@d){chomp $_;}
    return @d;
}

sub s{
    my $a=shift;#string to be split
    my $b=shift; #saparator
                 #default is \t
    if (!$b){
        $b="\t" ;
    }
    my @c=split(/$b/,$a);
    return @c;
}

