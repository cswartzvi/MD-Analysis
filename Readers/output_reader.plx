#! /usr/bin/perl -w

if( $#ARGV != 1 ){
   print "\n\nIncorrect number of arguments!\n\n";
   print " Argument 1) output file root (no number at the end).\n";
   print " Argument 2) Total number of output files to be considered\n\n\n";
   exit;
}
   
print " $ARGV[0] $ARGV[1]\n";


$output_root="$ARGV[0]";
$ncount = 0;

foreach $i (1 .. $ARGV[1]) {

   $output_file=$output_root.$i;
   print "$output_file \n";

   open INFILE, "<$output_file" ;
   @lines = <INFILE> ;
   close INFILE ;

   foreach $j ( 0 .. $#lines) {


      if ($lines[$j] =~ /Species  Temp/){
         $ncount += 1;

         @tmp1 = split /\s+/, $lines[$j + 1];
         @tmp2 = split /\s+/, $lines[$j + 2];
      
         $temp1[$ncount]=$tmp1[2]; 
         $temp2[$ncount]=$tmp2[2];
         $msd1[$ncount]=$tmp1[3];
         $msd2[$ncount]=$tmp2[3];

      }
   }
}  


   $temp_file="output_temps.dat";
   open OUTFILE, ">$temp_file";
   select OUTFILE ;

      foreach $i ( 1.. $ncount ) {

         printf("%6d  %6.2f  %6.2f \n", $i, $temp1[$i],  $temp2[$i] );

      }

   close OUTFILE;
   select STDOUT;


   $msd_file="output_msd.dat";
   open OUTFILE, ">$msd_file";
   select OUTFILE ;

      foreach $i ( 1.. $ncount ) {

         printf("%6d  %6.2f  %6.2f \n", $i, $msd1[$i],  $msd2[$i] );

      }

   close OUTFILE;
   select STDOUT;
