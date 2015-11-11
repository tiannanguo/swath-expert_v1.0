#!/usr/bin/perl
#tiannan,2015,IMSB,ETH
use Cwd;
use Getopt::Long;
use File::Slurp;
use threads;
use Storable;
use Data::Dumper;
use Statistics::LineFit;
use Statistics::Distributions;
$Data::Dumper::Sortkeys = sub { [sort {$a cmp $b} keys %{$_[0 ]}] };
use Array::transpose;

$in="nci60sw.txt";
$o="rt_norm_step2.sh";
@d=&oF($in);

open(OUT,">$o");
print OUT "mkdir rt_norm_com\n";
foreach my $d(@d){
          my @dd=&s($d);
          print OUT "bsub \"FileMerger -in rt_norm/$dd[1]\_*\.rtnorm.chrom.mzML -out rt_norm_com/$dd[1].rtnorm.chrom.mzML\"\n";

}
close OUT;




sub tg_anova3{
#One-way completely randomized
#follow tutorial at
#https://explorable.com/anova
#validated using online calculator
#http://turner.faculty.swau.edu/mathematics/math241/materials/anova/

#   use :
#   @d1=(2,3,7,2,6);
#   @d2=(2,3,7,5,10);
#   @d3=(3,2,4,3,5);
#   ($Fval,$f3,$f1,$pVal)=&tg_anova3(\@d1,\@d2,\@d3);
#   print "f is ",$Fval,"\np is ",$pVal,"\n";

    my $r=shift;
    my @d1=@$r;
    $r=shift;
    my @d2=@$r;
    $r=shift;
    my @d3=@$r;

    my $SS1_1=&tg_anova3_calSS1(@d1);  #sum of squares within groups
    my $SS1_2=&tg_anova3_calSS1(@d2);
    my $SS1_3=&tg_anova3_calSS1(@d3);

    my $SS1=$SS1_1+$SS1_2+$SS1_3; #within group squares

    my $SS2=&tg_anova3_calSS2(\@d1,\@d2,\@d3); #total squares

    my $SS3=$SS2-$SS1; #between group squares

    my $totalSampleNum=$#d1+$#d2+$#d3+3;
    my $groupNum=3;
    my $freedom1=$totalSampleNum-$groupNum;
    my $freedom3=$groupNum-1;

    my $Fval=($SS3/$freedom3)/($SS1/$freedom1);  #between group/within group

    #my $pVal=Statistics::Distributions::fdistr ($f3,f1,0.01);
    $pVal=Statistics::Distributions::fprob ($freedom3,$freedom1,$Fval);

    return ($Fval,$freedom3,$freedom1,$pVal);
}

sub tg_anova3_calSS1{
    my @d=@_;
    my $ave=&ave(@d);
    my $v=0;
    foreach (@d){
              $v+=($_-$ave)**2;
    }
    return $v;
}

sub tg_anova3_calSS2{
    my $r=shift;
    my @d1=@$r;
    $r=shift;
    my @d2=@$r;
    $r=shift;
    my @d3=@$r;
    my @d;
    foreach (@d1){push @d,$_;}
    foreach (@d2){push @d,$_;}
    foreach (@d3){push @d,$_;}
    my $ave=&ave(@d);
    my $v=0;
    foreach (@d){
              $v+=($_-$ave)**2;
    }
    return $v;
}


# randomly permutate @array in place
#fisher_yates_shuffle( \@array );    # permutes @array in place
sub fisher_yates_shuffle
{
    my $array = shift;
    my $i = @$array;
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}




sub digest_prot  {  #input protein sequence and cleavage

    my $prot=shift;

    $prot=~s/K/K\t/g;
    $prot=~s/R/R\t/g;
    my @pep=&s($prot);

    my @pep2=@pep;
    foreach (1..$#pep2){
             push @pep2,$pep[$_-1].$pep[$_];
    }
    foreach (2..$#pep2){
             push @pep2,$pep[$_-2].$pep[$_-1].$pep[$_];
    }

    @pep2=&unique(@pep2);

    my @pep3;
    foreach (@pep2){
              if (length($_)>=7 && length($_)<=55){
                  push @pep3,$_;
              }
    }
    @pep3=sort @pep3;

    return join(",",@pep3);  #return peptides array

}



sub get_6_element{
    my @d=@_;
    my @d2;
    foreach (1..6){
             $d2[$_-1]=$d[$_- 1];
    }
    foreach (1..6){
             shift @d;
    }
    return (\@d,\@d2);
}


sub lineFit{
    #use Statistics::LineFit;
    #more see http://search.cpan.org/dist/Statistics-LineFit-0.07/lib/Statistics/LineFit.pm
    my $r1=shift; #\@x
    my $r2=shift; #\@y
    my $lineFit=Statistics::LineFit->new();
    $lineFit->setData($r1,$r2);
    my ($a,$b)=$lineFit->coefficients();   #y=a+b*x
    my $rsq=$lineFit->rSquared();
    return ($a,$b,$rsq);
}

sub randomTimeNum{
    @timeData = localtime(time);
    my $timeNum=join( '', @timeData);
    my $randomNum=rand( 100000000);
    $randomNum=sprintf("%d" ,$randomNum);
    return $randomNum; #my $rtExtractFilename="mzML_rtExtract\_$timeNum\_$randomNum.txt";
}

sub getFiles{
    #eg: my @f=&getFiles("_cv.txt");
    my $namePattern=shift;
    $namePattern=~/\.([^\.]+)$/;
    my $end=$1;
    my @files=glob("*.$end");
    my @files2;
    foreach (@files){
              if (/$namePattern/){
                    push @files2,$_;
              }
    }
    return @files2;
}

sub median{
    @_ == 1 or die ( 'Sub usage: $median = median(\@array);' );
    my ($array_ref) = @_;
    my $count = scalar @$array_ref;
    # Sort a COPY of the array, leaving the original untouched
    my @array = sort { $a <=> $b } @$array_ref;
    if ($count % 2) {
    return $array[int($count/ 2)];
    } else {
    return ($array[$count/ 2] + $array[$count/ 2 - 1]) / 2;
    }
}


sub max{
    my @d=@_;
    my $max=shift @d;
    foreach (@d){
              if ($_>$max){
                   $max=$_;
              }
    }
    return $max;
}

sub min{
    my @d=@_;
    my $min=shift @d;
    foreach (@d){
              if ($_<$min){
                   $min=$_;
              }
    }
    return $min;
}

sub mean {
    my ($x)=@_;
    my $num = scalar(@{$x}) - 1;
    my $sum_x = '0';
    my $sum_y = '0';
    for (my $i = 1 ; $i < scalar(@{$x}); ++$i){
    $sum_x += $x->[$i][1 ];
    $sum_y += $x->[$i][2 ];
    }
    my $mu_x = $sum_x / $num;
    my $mu_y = $sum_y / $num;
    return($mu_x,$mu_y);
}
### ss = sum of squared deviations to the mean
sub ss {
    my ($x,$mean_x,$mean_y,$one,$two)=@_;
    my $sum = '0';
    for (my $i=1 ;$i<scalar(@{$x});++$i){
    $sum += ($x->[$i][$one]-$mean_x)*($x->[$i][$two]-$mean_y);
    }
    return $sum;
}
sub correlation {
    my ($x) = @_;
    my ($mean_x,$mean_y) = mean($x);
    my $ssxx=ss($x,$mean_x,$mean_y, 1, 1);
    my $ssyy=ss($x,$mean_x,$mean_y, 2, 2);
    my $ssxy=ss($x,$mean_x,$mean_y, 1, 2);
    my $correl=correl($ssxx,$ssyy,$ssxy);
    my $xcorrel=sprintf( "%.4f",$correl);
    return($xcorrel);
}
sub correl {
    my($ssxx,$ssyy,$ssxy)=@_;
    my $sign=$ssxy/abs($ssxy);
    my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
    return $correl;
}


sub oF2D{
    #usage:
    #my @d=&oF2D($in);
    #print "$d[0][0]\n";
    #
#example
#     my @d1=&oF2D($in1);
#     my @d2=&oF2D($in2);
#
#     #title
#     print OUT "$d1[0][0]";
#     foreach my $j (1..$#{$d1[0]}) {
#             print OUT "\t$d1[0][$j]";
#     }
#     print OUT "\n";
#
#
#     for my $i (1..$#d1) {
#           print OUT "$d1[$i][0]";
#           for my $j (1..$#{$d1[$i]}) {
#               $d1[$i][$j]=$d1[$i][$j]/$d2[$i][$j]*100;
#               print OUT "\t$d1[$i][$j]";
#           }
#           print OUT "\n";
#     }

    my $f=shift;
    my @d=&oF($f);
    my @d2;
    foreach my $d (@d){
            my @dd=&s($d);
            push @d2,\@dd;
    }
    return @d2;
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

sub oF{
    my $file=shift;
    my @d=read_file($file);
    foreach (@d){chomp $_;}
    return @d;
}

sub oF2{
    my $file=shift;
    open (IN,$file) || die "Error: can not open $file\n" ;
    my @d=<IN>;
    close IN;
    foreach (@d){chomp $_;}
    return @d;
}

sub unique{
    my @a=@_;
    @a=grep(($Last eq $_ ? 0 : ( $Last=$_, 1)),sort @a);
    return @a;
}

sub unique2{
    my @a=@_;
    my %a;
    foreach (@a){
              $a{$_}=1;
    }
    @a=keys %a;
    return @a;
}

sub str2array{
    my $p=shift;
    my @p;
    for (my $i=0 ;$i<length($p);$i++){
         push @p,substr($p,$i,1);
    }
    return @p;
}

sub printArray{
    my $r=shift;
    my @a=@$r;
    my $o=shift;
    open(OUT_printArray,">$o" );
    foreach (@a){
              print OUT_printArray "$_\n";
    }
    close OUT_printArray;
}


sub printHash{
    my $r=shift;
    my %a=%$r;
    my $o=shift;
    open(OUT_printHash,">$o" );
    foreach (sort {$a<=>$b} keys %a){
              print OUT_printHash "$_\t$a{$_}\n";
    }
    close OUT_printHash;
}



sub waitHr {
    my $hr=shift;
    $hr=$hr*3600 ;
    for (my $s=$hr;$s>0 ;$s--){
         print "$s of $hr\n" ;
         sleep(1 );
    }
}

sub ave{
    my @i=@_;
    my $a=0;
    foreach (@i){$a+=$_;}
    my $n=$#i+1;
    if ($n>0){
         $a/=$n;
         $a=sprintf("%.2f" ,$a);
    }
    else{
         $a=$i[0 ];
    }
    return $a;
}

sub sd{
    my @i=@_;
    my $ave=&ave(@i);
    my $sd;
    foreach (@i){
              $sd+=($_-$ave)**2;
    }
    $sd=($sd/($#i))**0.5 ;
    $sd=sprintf("%.2f" ,$sd);
    return $sd;
}

sub sum{
    my @i=@_;
    my $a=0;
    foreach (@i){$a+=$_;}
    return $a;
}


sub dos2unix{
    my $in=shift;
    my $o=shift;
    if (!$o){
        my $o= "tmp";
        open(INdos2unix,$in);
        my @d=<INdos2unix>;
        close INdos2unix;
        foreach (@d){
                  $_=~tr/\r\n//d;
        }
        open(OUTdos2unix,">$o" );
        foreach (@d){
                  print OUTdos2unix "$_\n";
        }
        close OUTdos2unix;
        unlink($in);
        rename ($o,$in);
    }
    else{
        open(INdos2unix,$in);
        my @d=<INdos2unix>;
        close INdos2unix;
        foreach (@d){
                  $_=~tr/\r\n//d;
        }
        open(OUTdos2unix,">$o" );
        foreach (@d){
                  print OUTdos2unix "$_\n";
        }
        close OUTdos2unix;
    }
}

sub divide_large_number{
    #modified 130614, use hash of hash
    #test
    #my %div=&divide_large_number(10,3);
    #foreach (sort keys %div){
    #         print "$_\t$div{$_}->{low}\t$div{$_}->{high}\n";
    #}
    #test this function using
    #$n=351;
    #$p=300;
    #%div=&read_file_lines($n,$p);
    #
    #foreach (sort {$a<=>$b} keys %div){
    #           my ($low,$high)=&s($div{$_});
    #           my $dif=$high-$low+1;
    #           print "$_\t$div{$_}\t$dif\n";
    #}
    my $n=shift;
    my $p=shift;
    my $res=$n%$p;
    my $n0=($n-$res)/$p;
    my $a=$res;      #print "a=$a\n";
    my $b=$n0+1;     #print "b=$b\n";
    my $c=$p-$res;   #print "c=$c\n";
    my $d=$n0;       #print "d=$d\n";

    my %div;
    my $i=1;
    foreach (1..$a){
             my %value;
             $value{low}=$i;
             $value{high}=$i+$b-1;
             $div{$_}=\%value;
             $i+=$b;
    }
    foreach ($a+1..$a+$c){
              my %value;
              $value{low}=$i;
              $value{high}=$i+$d-1;
              $div{$_}=\%value;
              $i+=$d;
    }
    return %div;
}


sub showTime {
    my ($sec,$min,$hr,$dayOfMonth,$month,$yearOffset,$dayOfWeek,$dayOfYear,$daylightSavings)=localtime;
    $hr='0' .$hr if ($hr< 10);
    $min='0' .$min if ($min< 10);
    $sec='0' .$sec if ($sec< 10);
    $month++;
    $month='0' .$month if ($month< 10);
    my $year=1900+$yearOffset;
    my $theTime="$year\.$month\.$dayOfMonth $hr:$min:$sec";
    print "time:" ,$theTime,"\n";
    return $theTime;
}

sub t_test_p_small{
    #bug: when degree of freedom goes >2000, it gets error, can not take log of -2.3434
    use Statistics::TTest;
             #how to use
             #my @a;
             #my @b;
             #my $p=&t_test_p(\@a,\@b);
    my $r1=shift;
    my $r2=shift;
    my @r1=@$r1;
    my @r2=@$r2;
    my $p="-";
    if (($#r1>1) and ($#r2>1 )){
        my $ttest = new Statistics::TTest;
        $ttest->set_significance(95);
        $ttest->load_data(\@r1,\@r2);
        $p=$ttest->{t_prob};
    }
    return $p;
}

sub t_test_p_paired{
    #this works for large set of values, but it only works for paried data
    use Statistics::DependantTTest;
    use Statistics::Distributions;
             #how to use
             #my @a;
             #my @b;
             #my $p=&t_test_p(\@a,\@b);
    my $r1=shift;
    my $r2=shift;
    my @x=@$r1;
    my @y=@$r2;
    my $p="-";
    if (($#x>1) and ($#y>1 )){
        my $ttest = new Statistics::DependantTTest;
        $ttest->load_data('x' ,@x);
        $ttest->load_data('y' ,@y);
        my ($tv,$df)=$ttest->perform_t_test( 'x','y' );    #print "tv $tv\ndf $df\n";
        $p=Statistics::Distributions::tprob($df,$tv);    #print "p $p\n";
        #convert one tail p value into 2 tail p value
        if ($tv>= 0){
             $p=$p*2;
        }
        else{
             $p=(1-$p)*2;
        }
    }
    return $p;
}