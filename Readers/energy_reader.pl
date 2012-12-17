#! /usr/bin/perl -w

   foreach $i (0 .. $#ARGV){
   
      if(lc($ARGV[$i]) eq "-df"){
         $data_file = $ARGV[++$i];
      }
      if(lc($ARGV[$i]) eq "-of"){
         $out_file = $ARGV[++$i];
      }
      if(lc($ARGV[$i]) eq "-tf"){
            $tot_frames = $ARGV[++$i];
      } 
      if(lc($ARGV[$i]) eq "-nat"){
         $nat = $ARGV[++$i];
      }  
   }
   
   open INFILE, "<$data_file" ;
   @lines = <INFILE> ;
   close INFILE ;

   printf("\n\n  Data File:   %s", $data_file);
   printf("\n  Total Frames: %f \n\n\n", $tot_frames);
   #printf("\n  Time Steps:   %f:\n\n", $dt);

   open OUTFILE, ">$out_file";
   select OUTFILE ;

   foreach $ns (0 .. ($tot_frames-1)){
      @temp =  split/\s+/, $lines[$ns];
      $KE   =  $temp[2];
      $ETOT =  $temp[5];
      $PE   =  $ETOT - $KE;

      #$t    =  $ns*$dt;

      printf("%8d %10.5f %10.5f %10.5f %10.5f %10.5f\n", $ns, $KE, $PE, $ETOT, ($PE/$nat), ($ETOT/$nat));

   }


   close OUTFILE;
   select STDOUT;
