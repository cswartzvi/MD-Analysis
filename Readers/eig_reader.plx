#! /usr/bin/perl -w

use POSIX;

#Eignvalue reader
# if PW then output file is the Main file
# if CP then output file id the *pos file

$eig_file=$ARGV[0];
$cal_type=$ARGV[1];
$states=$ARGV[2];

if ($#ARGV != 2 ) {
   print("\n  ERROR: Three needs to be 3 arugments\n");
   print(  "   1)Eigenvalue file\n");
   print(  "   2)Calculation type\n");
   print(  "   3)States \n\n");
   exit;
}

chomp($cal_type);

if($cal_type eq 'PW'){
   $ncol=8; 
   print ("\n  Calculations type: PW (must be the output file)\n");
} 
else{
   $ncol=10;
   print ("\n  Calculations type: CP (must be the *pos file) \n");
}

#Set rows, cols, and remains for reading eignevalues
$nrow  = ceil($states/$ncol);
$remain = $ncol - ($nrow*$ncol - $states);

open INFILE, "<$eig_file" ;
@eig_lines = <INFILE>;
close INFILE;

foreach $i ( 0 .. $#eig_lines) {


   #CP Calculation
   if($cal_type eq 'CP'){
      if ($eig_lines[$i] =~ /Eigenvalues/) {
         foreach $j ( 0 .. ($nrow-2) ) {
            @temp = split /\s+\n?/, $eig_lines[$i + 1 + $j];
            foreach $k ( 0 .. ($ncol-1) ) {
               $eigen[ $j*$ncol + $k ] = $temp[$k+1];
            }
         }
         @temp = split /\s+\n?/, $eig_lines[($i+1) + 2 + $nrow-2];
         foreach $k (0 .. ($remain-1) ) {
            $eigen[ ($nrow-1)*$ncol + $k ] = $temp[$k+1];
         }
      }
   }



   #PW Calculation
   if ($cal_type eq 'PW'){
      if ($eig_lines[$i] =~ /End of self-consistent calculation/) {
         foreach $j ( 0 .. ($nrow-2) ) {
            @temp = split /\s+\n?/, $eig_lines[$i + 4 + $j];
            foreach $k ( 0 .. ($ncol-1) ) {
               $eigen[ $j*$ncol + $k ] = $temp[$k+1];
            }
         }
         @temp = split /\s+\n?/, $eig_lines[($i+1) + 4 + $nrow-2];
         foreach $k (0 .. ($remain-1) ) {
            $eigen[ ($nrow-1)*$ncol + $k ] = $temp[$k+1];
         }
      }
   }


}


open FNAME, ">eig.dat";
select FNAME;

foreach $i (0 .. ($states-1)) {

   print "$eigen[$i] \n";
}
print "\n\n\n";

close FNAME;
select STDOUT;
