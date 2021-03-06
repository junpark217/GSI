#!/usr/bin/perl
#
# PROGRAM:get_gmi.pl  -
#  (1) acquires the GMI 1B and 1CR raw data  ( about  ~ 288  per day each )
#  (2) run gmi1cr_bufr.x  
#  (3)  archives the BUFR  files ( 1 per synoptic time)
#  (4)  "tar"  the raw GMI files for one day "tar" files( 1B and 1CR)
#  ( current day and pervious day)  and archive them.

# From Jianjun Jin Wargan  - gmi1cr_bufr.x 
#  Utilities to write GMI(1CR)  data into BUFR files

#

#
#
# 14 Octover   2014  Y. Kondratyeva

# The setting of the options and the module lookup paths will
# be done first using the BEGIN subroutine.  This section of the
# program executes before the rest of the program is even compiled.
# This way, a new path set via the -P option can be used to locate
# the modules to include at compile time while the remainder of the
# program is compiled.

# 14 December 2015  Y. Kondratyeva
#  If 'GMI_ACQUIRE_MACH' = MISS in Prep_Config.GET-GMI file , then we cp Input data 
# from Karki's directory.


# 4 April   2016  Y. Kondratyeva
#  New Option  "-g" added.
# if  defined $opt_g -->  $passive =  0  and transferred into Remote_utils.pm
#    in order  go ahead if No Input  data found. 


BEGIN {

# Keep track of errors within BEGIN block.

   $die_away = 0;
# Initialize output listing location

   $opt_O = 0;
# make env vars readily available
#--------------------------------
use Env qw( FORT_CONVERT20 );


# This module contains the getopts() subroutine.

   use Getopt::Std;
   use Getopt::Long;

# Get options and arguments

#  getopts('e:E:P:R:O:L:t:fab');

   GetOptions ( 'e=s',\$opt_e,
                'E=s',\$opt_E,
                'P=s',\$opt_P, 
                'R=s',\$opt_R,
                'O=s',\$opt_O,
                'L=s',\$opt_L,
                't=s',\$opt_t,
                'f',\$opt_f,
                'a',\$opt_a,
                'b',\$opt_b,
                'g',\$opt_g,
                'sched_cnfg=s',\$sched_cnfg,
                'sched_id=s',\$sched_id,
                'sched_synp=s',\$sched_synp,
                'sched_c_dt=s',\$sched_c_dt,
                'sched_dir=s',\$sched_dir,
                'sched_sts_fl=s',\$sched_sts_fl,
                'sched_hs=s',\$sched_hs );
#If get_gmi is initiated by the scheduler, construct table
# info. for "task_state" table of scheduler

   if ( defined( $sched_id ) )
   {
      $tab_status = 1;
      $tab_argv = "$sched_cnfg, $sched_id, $sched_synp, $sched_c_dt";
      $fl_name = "get_gmi";
      $comd_wrt = "$sched_dir/utility/status";
      $args = "$fl_name COMPLETE $tab_argv $sched_dir";
   }



# Processing environment

   if ( defined( $opt_e ) ) {
       $env = $opt_e;
   } else {
       $env = "ops";
   }

# The pre-processing configuration file.

   if ( defined( $opt_E ) ) {
      $PREP_CONFIG_FILE = $opt_E;
   } else {
      $PREP_CONFIG_FILE = "DEFAULT";
   }

  print "INPUT PREP_CONFIG_FILE = $PREP_CONFIG_FILE\n";

# Lag time for real-time processing (for llk mode only)

   if ( defined( $opt_L ) ) { 
      $LAG_TIME = $opt_L;
   } else {
      $LAG_TIME = 3;
   }

# Path to directory containing other GEOS DAS programs.
# Directory $GEOSDAS_PATH/bin will be searched for these
# programs.

   if ( defined( $opt_P ) ) { 
      $GEOSDAS_PATH = $opt_P;
   } else {
      $GEOSDAS_PATH = "DEFAULT";
   }

# Location of run-time configuration file.

   if ( defined( $opt_R ) ) { 
      $RUN_CONFIG_FILE = $opt_R;
   } else {
      $RUN_CONFIG_FILE = "DEFAULT";
   }

# Location of run-time configuration file.

   if ( defined( $opt_t ) ) {
      $syntime = $opt_t;
   } else {
      $syntime = "hd";
   }

           if  ($syntime eq '0') {
                $syntime = '00' ;
           }


          if   ($syntime eq '6') {
                $syntime = '06' ;
           }



# Processing environment

   if ( defined( $opt_g ) ) {
       $passive =  0 ;
   } else {
        $passive =  1 ;
   }




# ID for the preprocessing run.
#   We only run for 'flk'


$prep_ID=flk;


# print "INPUT prep_ID = $prep_ID\n";


# Location for output listings

   if ( $opt_O ) { 
      system ("mkdir -p $opt_O");
      if ( -w "$opt_O" ) {
        $listing_file    = "$opt_O/gmi_${prep_ID}.$$.listing";
        $listing_file_gz    = "$opt_O/gmi_${prep_ID}.$$.listing_gz";

        print "Standard output redirecting to $listing_file\n";
        open (STDOUT, ">$listing_file");
        open (STDERR, ">&" . STDOUT);
      }
      else {
        print "$0: WARNING: $opt_O is not writable for listing.\n";
      }
   }else{
        $listing_file = "STDOUT"
   }


# Set usage flag.  it is expected that 
# llk will have no synoptic time in it's argument list.
   
   $u_flag = 0;



 if ( $#ARGV < 0 || $ARGV[0] eq 'help' ) {

          $u_flag = 1;
      }

   if ( $u_flag == 1 ) {
       print STDERR <<'ENDOFHELP';
Usage:

   get_gmi.pl [-e] [-E Prep_Config] [-P GEOSDAS_Path] [-R Run_Config] [ -O output_location ] [-L lag_time ] [-t synoptic time] [ process_date ]

   Normal options and arguments:

   -e    Processing environment (default = ops)

   -f    Force flag.If no NETCDF data for a synoptic time, issue error,
         but keep on processing.If no data on remote machine(rget) for
         a synoptic time,or No input data for converter, - issue error,
         but keep on processing.

   -d   Force flag.  If "wget" exit with error ,or no files(or size=0)
        for yesterday or today day on remote machine(wget)
         --> continue processing.

   -a    ARCHIVE flag.If there is option "a",then we GZIP and ARCHIVE INPUT ORBIT FILES 

   -b    REBLOCK flag.If there is option "b",then we REBLOCK BUFR files before archive.


   -O output_location
         This is the full path to the output listings (both STDERR and STDOUT).

   -E Prep_Config
         The full path to the preprocessing configuration file.  This file contains
         parameters needed by the preprocessing control programs. If not given, a
         file named $HOME/$prep_ID/Prep_Config is used.  get_gmi.pl exits with an
         error if neither of these files exist.

         The parameters set from this file are
   -t synoptic time
         This is the synoptic time to process.


         WORK_DIR
            The top level working directory in which to run the preprocessing.  A
            subdirectory with the name of the preprocessing ID is made in this 
            directory (i.e., $WORK_DIR/$prep_ID/$process_date), and the work is done 
            there.


         GMI_STAGE_DIR
            The location in which to stage the BFR files for use by the DAS.

         GMI_ARCHIVE_LOC
            The location in which to archive the BUFR and GMI input  files.  Files will be 
            stored in subdirectories according to their Sat_ID  and valid date.  As an 

         GMI_ACQUIRE_MACH
            This is ftp location for the raw native GMI
            data.


         GMI_ACQUIRE_PATH_1B
            This is the file path template for the raw native GMI 1B
            data.

         GMI_ACQUIRE_PATH_1CR
            This is the file path template for the raw native GMI 1CR
            data.


         GMI_BUFR
            This is the output BUFR filename template.


         GMI_TABLE_DIR
            Directory where bufr tables and other resource files used by the
            processing are stored.

         Satellite identification tag for this run of the GMI pre-processing.

   prep_ID
         Identification tag for this run of the GMI pre-processing.

   process_date
         Date in YYYYMMDD format to process.  If not given, then today's date (in
         terms of GMT) will be processed.

   synoptic time 
         Synoptic time  in  hh format to process (00,06,12,18 or 0,6,12,18).  If not given, then  the whole day is processed ( $process date)

   Options useful in development mode.  These are not necessary (and should not be
   used) when running this program in the usual operational environment.

   -P GEOSDAS_Path
         Path to directory containing other GEOS DAS programs.  The path is 
         $GEOSDAS_PATH, where $GEOSDAS_PATH/bin is the directory containing these
         programs.  If -P GEOSDAS_Path is given, then other required programs not 
         found in the directory where this program resides will be obtained from 
         subdirectories in GEOSDAS_Path - the subdirectory structure is assumed 
         to be the same as the operational subdirectory structure.  The default is 
         to use the path to the subdirectory containing this program, which is what 
         should be used in the operational environment.

   -R Run_Config
         Name of file (with path name, if necessary) to read to obtain the 
         run-time (execution) configuration parameters.  get_gmi.pl reads this
         file to obtain configuration information at run time.  

         If given, get_gmi.pl uses this file.  Otherwise, get_gmi.pl looks for a 
         file named "Run_Config" in the user's home directory, then the 
         $GEOSDAS_PATH/bin directory.  $GEOSDAS_PATH is given by the -P option if
         set, or it is taken to be the parent directory of the directory in which this
         script resides.  It is an error if get_gmi.pl does not find this file, 
         but in the GEOS DAS production environment, a default Run_Config file is always 
         present in the bin directory.

    -L   Lag Time
         This option is to be used in llk real time mode only.  This is the lag time, in 
         days, before the today date.

ENDOFHELP
      $die_away = 1;
   }


# This module locates the full path name to the location of this file.  Variable
# $FindBin::Bin will contain that value.

   use FindBin;

# This module contains the dirname() subroutine.

   use File::Basename;

# If default GEOS DAS path, set path to parent directory of directory where this
# script resides.  

   if ( $GEOSDAS_PATH eq "DEFAULT" ) {
      $GEOSDAS_PATH = dirname( $FindBin::Bin );
   }

# Set name of the bin directory to search for other programs needed by this one.

   $BIN_DIR = "$GEOSDAS_PATH/bin";

# Get the name of the directory where this script resides.  If it is different 
# than BIN_DIR, then this directory will also be included in the list of 
# directories to search for modules and programs.

   $PROGRAM_PATH = $FindBin::Bin;

# Now allow use of any modules in the bin directory, and (if not the same) the
# directory where this program resides.  (The search order is set so that
# the program's directory is searched first, then the bin directory.)

   if ( $PROGRAM_PATH ne $BIN_DIR ) {
      @SEARCH_PATH = ( $PROGRAM_PATH, $BIN_DIR );
   } else {
      @SEARCH_PATH = ( $BIN_DIR );
   }

# Set module environment for Fortran executable.

   print "source g5_modules.\n";
   do "${BIN_DIR}/g5_modules_perl_wrapper";

}	# End BEGIN

# Any reason to exit found during the BEGIN block?

if ( $die_away == 1 ) {
   exit 1;
}

# Include the directories to be searched for required modules.

use lib ( @SEARCH_PATH );

# Set the path to be searched for required programs.

$ENV{'PATH'} = join( ':', @SEARCH_PATH, $ENV{'PATH'} );

# This module contains the extract_config() subroutine.
use Extract_config;

# Archive utilities: gen_archive
use Arch_utils;

# This module contains the z_time(), dec_time() and date8() subroutines.
use Manipulate_time;

# Error logging utilities.
use Err_Log;

# Record FAILED to schedule status file.
use Recd_State;

# This module contains the mkpath() subroutine.

use File::Path;
use File::Copy;

# This module contains the rget() routine.

use Remote_utils;

# This module contains the julian_day subroutine.

use Time::JulianDay;

#Initialize exit status

$exit_stat = 0;

# Set Event/Error log message prefix.

    if ( defined( $sched_id ) )  {
       $err_pref="$sched_id";
    }
    elsif ( ( ${opt_n} ) ) {
       $err_pref="get_gmi_${opt_n}";
    }
    else {
       $err_pref="get_gmi";
    }


# Write start message to Event Log

err_log (0, "get_gmi.pl", "$prep_ID","$env","-1",
        {'err_desc' => "${err_pref}: prep_ID get_gmi.pl job running for has started - Standard output redirecting to $listing_file"});

# Use Prep_Config file under the preprocessing run's directory in the user's home directory
# as the default.

if ( "$PREP_CONFIG_FILE" eq "DEFAULT" ) {
   $PREP_CONFIG_FILE = "$ENV{'HOME'}/$prep_ID/Prep_Config";
}

# Does the Prep_Config file exist?  If not, die.
if ( ! -e "$PREP_CONFIG_FILE" ) {
    err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
	     {'err_desc' => "${err_pref}: error $PREP_CONFIG_FILE not found while running for ."});
    recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
    die "error $PREP_CONFIG_FILE not found.";
}

# If date given, use that, 
# otherwise use today's date (GMT).


    if ( $#ARGV >= 0 ) {
	$process_date = date8( $ARGV[0] );
    } else {
	
# Get current date (YYYYMMDD) in GMT, and set the process date to be 
# $LAG_TIME days prior.  
	
	$process_date = ( z_time() )[0];
	($process_date, $process_time) = inc_time ($process_date, 0, -$LAG_TIME, 0);
    }

# The date strings in the error messages and on the listing files is a function of 
# the mode we're running.

    $err_time = "${process_date}";

   if ( $syntime eq 'hd') {
            $err_time = "${process_date}";
   } else {
           $err_time = "${process_date}.${syntime}z";
    }





# Find the locations in which to stage and archive the BUFR files.
( $GMI_STAGE_DIR = extract_config( "GMI_STAGE_DIR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_STAGE_DIR configuration value\n";

# Get the location  to archive GMI data( BUFR and INPUT) 

( $GMI_ARCHIVE_LOC = extract_config( "GMI_ARCHIVE_LOC", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_ARCHIVE_LOC configuration value\n";


# Get the location, directory, and file names for the Input  GMI data.

 ( $GMI_ACQUIRE_MACH = extract_config( "GMI_ACQUIRE_MACH", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
    or die "(get_gmi.pl) ERROR - can not set GMI_ACQUIRE_MACH configuration value\n";


( $GMI_ACQUIRE_PATH_1B  = extract_config( "GMI_ACQUIRE_PATH_1B", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_ACQUIRE_PATH_1B  configuration value\n";

$template_path_1B=$GMI_ACQUIRE_PATH_1B ;

( $GMI_ACQUIRE_PATH_1CR  = extract_config( "GMI_ACQUIRE_PATH_1CR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_ACQUIRE_PATH_1CR  configuration value\n";

$template_path_1CR=$GMI_ACQUIRE_PATH_1CR ;


# Get the name of the working directory for the observation preprocessing.

( $GMI_WORK_DIR = extract_config( "GMI_WORK_DIR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_WORK_DIR configuration value\n";


#   Get the template name of TABLES

   ( $GMI_TABLE_DIR = extract_config( "GMI_TABLE_DIR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
      or die "(get_gmi.pl) ERROR - can not set GMI_TABLE_DIR configuration value\n";


# Get the template for the output GMI bufr file.

( $GMI_BUFR = extract_config( "GMI_BUFR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
   or die "(get_gmi.pl) ERROR - can not set GMI_BUFR configuration value\n";

$template_BUFR=$GMI_BUFR ;


# --------------------------------------
# Assign - assigns file name to Fortran units.

sub Assign {

    my ( $fname, $lu ) = @_;

    $f77name = "fort.$lu";
    unlink($f77name) if ( -e $f77name ) ;
    symlink("$fname","$f77name");

  }
# -------------------------------------


#!/usr/bin/perl

#-- SUBROUTINE GET GMI BY SYNOPTIC ---------------------
#   Subroutine go through list of GMI files and sort them by synoptic time
#   Write lists of files  by synoptic time .
#   Create array from names of files by synoptic time
#   Return list of files,arrays of file's names and count of files by synoptic time


sub get_bysyn_gmi {
   my (  $date_yesterday,$date_today ) = @_;

#  $date_yesterday,$date_today  are in FORM YYYYMMDD  
            $k00=  0;
            $k06=  0;
            $k12=  0;
            $k18=  0;
    $list00 =' ';
    $list06 =' ';
    $list12 =' ';
    $list18 =' ';


    @lst_file00='';
    @lst_file06='';
    @lst_file12='';
    @lst_file18='';
      

#  Files for GMI 1B can be like

#   0123456789012345678901234567890123456789012345678901
#   1B.GPM.GMI.TB2014.20141006-S195647-E200145.V03B.RT-H5

#  Files for GMI 1CR  can be like

#   0123456789012345678901234567890123456789012345678901
#   1C-R.GPM.GMI.XCAL2014-N.20141006-S200147-E200645.V03B.RT-H5
#  ????????????????????//

# FOR GMI  1B
#   0123456789012345678901234567890123456789012345678901
#   1B.GPM.GMI.TB2014.20141006-S195647-E200145.V03B.RT-H5

#     $date_current= =substr( $gos,18,8) ;
#     $gap_stime=substr( $gos,28,2) ;
#     $gap_etime=substr( $gos,36,2) ;

# FOR GMI  1C-R
#   0123456789012345678901234567890123456789012345678901
#   1C-R.GPM.GMI.XCAL2014-N.20141006-S200147-E200645.V03B.RT-H5

#     $date_current= =substr( $gos,24,8) ;
#     $gap_stime=substr( $gos,34,2) ;
#     $gap_etime=substr( $gos,42,2) ;



      while ( defined($nextname =  <1C-R.GPM.GMI.XCAL*>)) {

     $nextname =~s#.*/##;  # remove part before last  slash

      $gos ="$nextname";

#===================================
# FOR GMI  1C-R
#   0123456789012345678901234567890123456789012345678901
#   1C-R.GPM.GMI.XCAL2014-N.20141006-S200147-E200645.V03B.RT-H5

      $date_current =substr( $gos,24,8) ;
      $gap_time=substr( $gos,34,3) ;



  if ( $date_current == $date_yesterday) {

#   S200147   ==>   $gap_time = 200
#   S030147   ==>   $gap_time = 030
#    0z   Yesterday   Current time >= 205
#         TODAY Current    time < 30 ( time < 030)

      if  ($gap_time >= 205  ) {

         $list00='$list00 $nextname';
         $lst_file00 [$k00] = $nextname ;
         $k00= $k00 +1;
   
       }
    
# end of date_yesterday
     }

 if ( $date_current == $date_today) {

       if ( $gap_time < 30) {
         $list00='$list00 $nextname';
         $lst_file00 [$k00] = $nextname ;
         $k00= $k00 +1;
        }


#    6z   Current  time >= 25 ( time = 025)
#        Current    time < 90 ( time < 090)

    if ( $gap_time < 90 && $gap_time >= 25 ) {

         $list06="$list06 $nextname";
          $lst_file06 [$k06] = $nextname ;
         $k06= $k06 +1;
     }

#    12z   Current  time >= 085( time =085)
#        Current    time < 150 


     if ( $gap_time < 150 && $gap_time >= 85) {


         $list12="$list12 $nextname";
          $lst_file12 [$k12] = $nextname ;
         $k12= $k12 +1;
     }

#    12z   Current  time >= 145
#        Current    time < 210 

     if ( $gap_time < 210 && $gap_time >= 145) {
         $list18="$list18 $nextname";
          $lst_file18 [$k18] = $nextname ;
         $k18= $k18 +1;
     }

# end of today
}

#  END  LOOP BY inlist ( input files)
   }


       return ($k00,$k06,$k12,$k18,$lst_file00,$lst_file06,$lst_file12,$lst_file18,$list00,$list06,$list12,$list18);

   }

#--END OF SUBROUTINE GET BY SYNOPTIC ---------------------
# ---------------------------------------

# =====================================================================
#   SUBROUTINE: CREATE LIST of input FILES  by day . Return lengh of list.
# =====================================================================
     sub get_native {

      my (${filename_beg}, $date) = @_;
              $real_len_full=0 ;

# Name is 1B.GPM.GMI.TB.%y4%m2%d2*H5   for 1B files
# Name is 1C-R.GPM.GMI.XCAL.%y4%m2%d2*H5     for 1CR files
#    while ( defined($nextname =  <1B.GPM.GMI.TB*$date*H5>)) 
# or
#     while ( defined($nextname =  <1C-R.GPM.GMI.XCAL*$date*H5>))


    while ( defined($nextname =  <${filename_beg}*${date}*H5>))  {

         $nextname =~s#.*/##;  # remove part before last  slash
         $full_list[$real_len_full] =$nextname;
#               print "VNUTRI full_list  = $full_list[$real_len_full] \n";


          $real_len_full=${real_len_full} +1;


#  END  LOOP BY  input files
    }
          print " real_len_full =  $real_len_full\n";

          return ($real_len_full,@full_list);

   }

# =========================================
#  END of  SUBROUTINE get_native
# =========================================


# ---------------------------------------

############################################
# MAKE WORK_DIR  
################################################################
# Get the work path and Make it.  (mkpath default mode is 0777, which is what we want.)


# Change into WORK directory it.  and clean it



#  **********************************************************

   $gmi_work="$GMI_WORK_DIR/$prep_ID/${process_date}";

if ( ! -d ${gmi_work} ) {

    unless (defined eval {mkpath( "${gmi_work}" )})
    {
       err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
               {'err_desc' => "${err_pref}: Fatal Error: Unable to make directory ${gmi_work}."});
       recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
       die "Cannot make ${gmi_work}";
    }
}


# Change into WORK directory it.

  chdir "${gmi_work}" or die "Cannot cd to ${gmi_work}: $!\n";


             $rc=system("rm -f * ");

#  **********************************************************
# Make  STAGE  directories if they don't already exist.


$GMI_STAGE_DIR = token_resolve( "$GMI_STAGE_DIR", $process_date, $synhour );

if ( ! -d "$GMI_STAGE_DIR" ) {
#   mkpath( "$STAGE_DIR" ) or die "Cannot make $STAGE_DIR";

    unless (defined eval {mkpath( "${GMI_STAGE_DIR}" )})
    {
       err_log (4, "get_gmi.pl", "$err_time","-1",
               {'err_desc' => "${err_pref}: Fatal Error: Unable to make directory ${GMI_STAGE_DIR}."});
       recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
       die "Cannot make ${GMI_STAGE_DIR}";
    }

}



################################################################


# -----------------------------------------------
#   Started to ftp input GMI 1B  files , run  
# -----------------------------------------------

# We need files from the current and the previous day to create  BUFR file 
# (4 synoptic times of GMI data). 
# First we'll get the list of today's file.


  ( $process_date_m1, $current_time ) = inc_time ($process_date, $current_time, -1, 0);

# -----------------------------------------------
#   Started to ftp input GMI 1B  files , run  
# -----------------------------------------------

# We need files from the current and the previous day to create  BUFR file 
# (4 synoptic times of GMI data). 
# First we'll get the list of today's file.

######################################################
#         GMI 1B 
#   Get files GMI 1B for TODAY
##################################################
   $GMI_ACQUIRE_PATH=$template_path_1B;
  $GMI_ACQUIRE_PATH =token_resolve("${GMI_ACQUIRE_PATH}",$process_date);
     print "POPALI INPUT GMI_ACQUIRE_PATH  = $GMI_ACQUIRE_PATH \n";


@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

# Paste the array back together into a string leaving off the filename.

${GMI_ACQUIRE_DIR} = join('/',@list[0 .. $list_len-2]);

# Get the list of files available on the remote server and then grab those files.
#- for current day   -  can be 288 files 

  $filename="${GMI_ACQUIRE_FILE}";
   $filename_today_1B = $filename ;
    print "POPALI INPUT LIST filename = $filename\n";
# pppppppppppppppppppppppppppppppppppppppppppppppppppppppp

# Check   - from where to take Input files.

     if ($GMI_ACQUIRE_MACH  eq  'MISS' ) {

#   Copy Input files from Karki's directory.
         $rc=system("cp $GMI_ACQUIRE_PATH   . ");

      }
              else  {

	#   Copy (rget) Input files from ftp site .

	# Get the list of files available on the remote server and then grab those files.
	#- for current day   -  can be 288 files 

#        $passive = 0;
      print "PERED VAJNO passive  = $passive \n";



# %options = $passive ;   # ne rabotaet

  %options = ('passive'  => $passive) ;
	  print "Calling rflist: $GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename\n";
# chomp(@remote_list_today=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename"));

  chomp(@remote_list_today=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename",\%options));


	# print " remote_list_today  = @remote_list_today  .\n";

	$remote_list_today_len=@remote_list_today;
	  print " LENGHT remote_list_today  = $remote_list_today_len  .\n";


	#Check to see if there is file to grab
	if ($remote_list_today_len > 0){

        
	$i=0;

           $need_file_1B = $remote_list_today[$i];

           print " NADONADO 1-i FILE  remote_list_today = $remote_list_today[i] \n";


#     print     $remote_list_today[$i]
	   while  ($i < $remote_list_today_len ){

	   $remote_namepath="$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$remote_list_today[$i]";


	$remote_list_today_len=@remote_list_today;

	    $rget_retcode = rget("$remote_namepath", "$remote_list_today[$i]",
		  {'debug' => "1",'run_config' => "$RUN_CONFIG_FILE" });

	    if( $rget_retcode != 1) {
		err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
			 {'err_desc' => "${err_pref}: Error: rget ftp failed trying to get $remote_list_today[$i] from  $remote_namepath."});
	     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
		die "Error: rget ftp failed trying to get $remote_list_today[$i] from $remote_namepath. ";
	     }

#--VOTVOT------------------------

# end if $i < $remote_list_today_len 
		$i=$i+1;
	    }

	}
	#  If NO FILE for today , $remote_list_today_len <=  0
	    else {

# if ($remote_list_today_len <= 0)
	       err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
	       {'err_desc' => "${err_pref}: NON-Fatal Warning: There are no files available on the remote machine for processing GMI 1B on $process_date."});
		$archive_err ++;

	   print "POPALI UJAS 10  remote_list_today_len = $remote_list_today_len\n";
	       if ( ! $opt_f ) {
#     print "STOP : NO OPTION F \n";
	    recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
		  die "Error:  There be no INPUT  data available for date: $process_date  .  This error occurred while processing GMI 1B data .";
	       }

	    }

#  END of    Check   - from where to take Input files.
    }
# pppppppppppppppppppppppppppppppppppp
################################################################


# -----------------------------------------------
#   Started to cp input GMI 1CR  files 
# -----------------------------------------------

# We need files from the current and the previous day to create  BUFR file 
# (4 synoptic times of GMI data). 
# First we'll get the list of today's file.

######################################################
#         GMI 1CR 
#   Get files GMI 1CR for TODAY
##################################################
   $GMI_ACQUIRE_PATH=$template_path_1CR;
  $GMI_ACQUIRE_PATH =token_resolve("${GMI_ACQUIRE_PATH}",$process_date);


@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

# Paste the array back together into a string leaving off the filename.

${GMI_ACQUIRE_DIR} = join('/',@list[0 .. $list_len-2]);

# Get the list of files available on the remote server and then grab those files.
#- for current day   -  can be 288 files 

  $filename="${GMI_ACQUIRE_FILE}";
     $filename_today_1CR = $filename ;

    print "POPALI INPUT LIST filename = $filename\n";

# pppppppppppppppppppppppppppppppppppppppppppppppppppppppp

# Check   - from where to take Input files.

       if ($GMI_ACQUIRE_MACH  eq  'MISS' ) {

#   Copy Input files from Karki's directory.
             $rc=system("cp $GMI_ACQUIRE_PATH   . ");

       }
          else {

#   Copy (rget) Input files from ftp site .

# Get the list of files available on the remote server and then grab those files.
#- for current day   -  can be 288 files

# %options = $passive ;
# %options = ('passive'  => 0) ;
      print "PERED VAJNO passive  = $passive \n";
  %options = ('passive'  => $passive) ;


  print "Calling rflist: $GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename\n";
# chomp(@remote_list_today=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename"));

  chomp(@remote_list_today=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename",\%options));


# print " remote_list_today  = @remote_list_today  .\n";

$remote_list_today_len=@remote_list_today;
  print " LENGHT remote_list_today  = $remote_list_today_len  .\n";


#Check to see if there is file to grab
if ($remote_list_today_len > 0){

$i=0;
           $need_file_1CR = $remote_list_today[$i];
   while  ($i < $remote_list_today_len ){

   $remote_namepath="$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$remote_list_today[$i]";


$remote_list_today_len=@remote_list_today;

    $rget_retcode = rget("$remote_namepath", "$remote_list_today[$i]",
          {'debug' => "1",'run_config' => "$RUN_CONFIG_FILE" });

    if( $rget_retcode != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: Error: rget ftp failed trying to get $remote_list_today[$i] from  $remote_namepath."});
     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
        die "Error: rget ftp failed trying to get $remote_list_today[$i] from $remote_namepath. ";
     }


# end if $i < $remote_list_today_len 
        $i=$i+1;
    }

}
#  If NO FILE for today , $remote_list_today_len <=  0
    else {

# if ($remote_list_today_len <= 0)
       err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
       {'err_desc' => "${err_pref}: NON-Fatal Warning: There are no files available on the remote machine for processing GMI 1CR on $process_date."});
        $archive_err ++;

   print "POPALI UJAS 10  remote_list_today_len = $remote_list_today_len\n";
       if ( ! $opt_f ) {
    recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error:  There be no INPUT  data available for date: $process_date  .  This error occurred while processing GMI 1CR data .";
       }

    }

#  END of    Check   - from where to take Input files.
# ?????
    }
# ==========================================

  ( $process_date_m1, $current_time ) = inc_time ($process_date, $current_time, -1, 0);

#  444444444444444444444444444444444444444444444444444444444444444
     if ( ($syntime eq '00') or ($syntime eq 'hd')) {
#  If we have to get 0z , or the whole day ( syntime='00',or='0', or syntime = 'hd')
# We need files from the current and the previous day to create  BUFR file 
# (4 synoptic times of GMI data). 


################################################################


# THEN  we'll get the list of yesterday's file.

######################################################
#         GMI 1B 
#   Get files GMI 1B for YESTERDAY
##################################################
( $process_date_m1, $current_time ) = inc_time ($process_date, $current_time, -1, 0);


   $GMI_ACQUIRE_PATH=$template_path_1B;
  $GMI_ACQUIRE_PATH =token_resolve("${GMI_ACQUIRE_PATH}",$process_date_m1);


@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

# Paste the array back together into a string leaving off the filename.

${GMI_ACQUIRE_DIR} = join('/',@list[0 .. $list_len-2]);

# Get the list of files available on the remote server and then grab those files.
#- for a day   -  can be 288 files 

  $filename="${GMI_ACQUIRE_FILE}";
  $filename_yesterday_1B = $filename ;
    print "POPALI INPUT LIST filename = $filename\n";

# pppppppppppppppppppppppppppppppppppppppppppppppppppppppp

# Check   - from where to take Input files.

      if ($GMI_ACQUIRE_MACH  eq  'MISS' ) {

#   Copy Input files from Karki's directory.

   $rc=system("cp $GMI_ACQUIRE_PATH   . ");

      }
          else {
#   Copy (rget) Input files from ftp site .

# Get the list of files available on the remote server and then grab those files.
#- for current day   -  can be 288 files

# %options = $passive ;
      print "PERED VAJNO passive 3 = $passive \n";

  %options = ('passive'  => $passive) ;


    print "POPALI GMI_ACQUIRE_DIR = $GMI_ACQUIRE_DIR \n";

  print "Calling rflist: $GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename\n";
# chomp(@remote_list_yesterday=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename"));

  chomp(@remote_list_yesterday=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename",\%options));


# print " remote_list_yesterday  = @remote_list_yesterday  .\n";

$remote_list_yesterday_len=@remote_list_yesterday;
  print " LENGHT remote_list_yesterday  = $remote_list_yesterday_len  .\n";


#Check to see if there is file to grab
if ($remote_list_yesterday_len > 0){

$i=0;
  $need_file_1B = $remote_list_yesterday[$i];
   while  ($i < $remote_list_yesterday_len ){

   $remote_namepath="$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$remote_list_yesterday[$i]";


$remote_list_yesterday_len=@remote_list_yesterday;

    $rget_retcode = rget("$remote_namepath", "$remote_list_yesterday[$i]",
          {'debug' => "1",'run_config' => "$RUN_CONFIG_FILE" });

    if( $rget_retcode != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: Error: rget ftp failed trying to get $remote_list_yesterday[$i] from  $remote_namepath."});
     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
        die "Error: rget ftp failed trying to get $remote_list_yesterday[$i] from $remote_namepath. ";
     }


# end if $i < $remote_list_yesterday_len 
        $i=$i+1;
    }

}
#  If NO FILE for yesterday , $remote_list_yesterday_len <=  0
    else {

# if ($remote_list_yesterday_len <= 0)
       err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
       {'err_desc' => "${err_pref}: NON-Fatal Warning: There are no files available on the remote machine for processing GMI 1B on $process_date_m1."});
        $archive_err ++;

   print "POPALI UJAS 10  remote_list_yesterday_len = $remote_list_yesterday_len\n";
       if ( ! $opt_f ) {
    recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error:  There be no INPUT  data available for date: $process_date_m1  .  This error occurred while processing GMI 1B data .";
       }

    }

#  END of    Check   - from where to take Input files.
# ?????
    }



#44444444444444444444444444444444444444444444444444444444444


#  444444444444444444444444444444444444444444444444444444444444444
#  If we have to get 0z , or the whole day ( syntime = '00', or syntime = 'hd')
# We need files from the current and the previous day to create  BUFR file 
# (4 synoptic times of GMI data). 


################################################################


# THEN  we'll get the list of yesterday's file.

######################################################
#         GMI 1CR 
#   Get files GMI 1CR for YESTERDAY
##################################################
   $GMI_ACQUIRE_PATH=$template_path_1CR;
  $GMI_ACQUIRE_PATH =token_resolve("${GMI_ACQUIRE_PATH}",$process_date_m1);


@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

# Paste the array back together into a string leaving off the filename.

${GMI_ACQUIRE_DIR} = join('/',@list[0 .. $list_len-2]);

# Get the list of files available on the remote server and then grab those files.
#- for a day   -  can be 288 files 

  $filename="${GMI_ACQUIRE_FILE}";
    $filename_yesterday_1CR = $filename ;
    print "POPALI INPUT LIST filename = $filename\n";
# pppppppppppppppppppppppppppppppppppppppppppppppppppppppp

# Check   - from where to take Input files.

      if ($GMI_ACQUIRE_MACH  eq  'MISS' ) {

#   Copy Input files from Karki's directory.

         $rc=system("cp $GMI_ACQUIRE_PATH   . ");

       }
            else {
#   Copy (rget) Input files from ftp site .

# Get the list of files available on the remote server and then grab those files.
#- for current day   -  can be 288 files

# %options = $passive ;

      print "PERED VAJNO passive  4 = $passive \n";



  %options = ('passive'  => $passive) ;


  print "Calling rflist: $GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename\n";
# chomp(@remote_list_yesterday=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename"));
  chomp(@remote_list_yesterday=rflist("$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$filename",\%options));

# print " remote_list_yesterday  = @remote_list_yesterday  .\n";

$remote_list_yesterday_len=@remote_list_yesterday;
  print " LENGHT remote_list_yesterday  1CR = $remote_list_yesterday_len  .\n";


#Check to see if there is file to grab
if ($remote_list_yesterday_len > 0){

$i=0;
  $need_file_1CR = $remote_list_yesterday[$i];
   while  ($i < $remote_list_yesterday_len ){

   $remote_namepath="$GMI_ACQUIRE_MACH:$GMI_ACQUIRE_DIR/$remote_list_yesterday[$i]";


$remote_list_yesterday_len=@remote_list_yesterday;
# print "POPALI UJAS 4  remote_list_yesterday_len = $remote_list_yesterday_len\n";

    $rget_retcode = rget("$remote_namepath", "$remote_list_yesterday[$i]",
          {'debug' => "1",'run_config' => "$RUN_CONFIG_FILE" });

    if( $rget_retcode != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: Error: rget ftp failed trying to get $remote_list_yesterday[$i] from  $remote_namepath."});
     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
        die "Error: rget ftp failed trying to get $remote_list_yesterday[$i] from $remote_namepath. ";
     }

# end if $i < $remote_list_yesterday_len 
        $i=$i+1;
    }

}
#  If NO FILE for yesterday , $remote_list_yesterday_len <=  0
    else {

# if ($remote_list_yesterday_len <= 0)
       err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
       {'err_desc' => "${err_pref}: NON-Fatal Warning: There are no files available on the remote machine for processing GMI 1CR on $process_date_m1."});
        $archive_err ++;

   print "POPALI UJAS 10  remote_list_yesterday_len = $remote_list_yesterday_len\n";
       if ( ! $opt_f ) {
#     print "STOP : NO OPTION F \n";
    recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error:  There be no INPUT  data available for date: $process_date_m1  .  This error occurred while processing GMI 1CR data .";
       }

    }
#  END of    Check   - from where to take Input files.
    }
# ppppppppppppppppppppppppppppppppppppppppppppp

#   End if we  need the yesterday files

  }
#44444444444444444444444444444444444444444444444444444444444

# =======================================================

#        $date_today = $process_date ;
#           $date_yesterday = $process_date_m1 ;

       get_bysyn_gmi ( $process_date_m1,$process_date);

#      print "KONEZ  m=    is $m\n";

       print "AFTER get_bysyn  k00=    is $k00\n";
       print "AFTER get_bysyn  k06=    is $k06\n";

       print "AFTER get_bysyn  k12=    is $k12\n";
       print "AFTER get_bysyn  k18=    is $k18\n";
 

#     print "AFTER get_bysyn  list00=    is $list00\n";
#    print "AFTER get_bysyn lst_file06  =  @lst_file06 \n";

#    print "AFTER get_bysyn lst_file06 (0) =  @lst_file06[$kk0]\n";
#    print "AFTER get_bysyn lst_file06 (1) =  @lst_file06[$kk1]\n";
#    print "AFTER get_bysyn lst_file06 (2) =  @lst_file06[$kk2]\n";
#   print "     \n";
#    print "AFTER get_bysyn lst_file06 (3) =  @lst_file06[$kk3]\n";



##########################
#  Copy BUFR TABLE for GMI
###########################
   $rc=system("cp $GMI_TABLE_DIR/GMI_bufr_table_1CR ${gmi_work}  ");
##########################



   if ( $syntime eq 'hd') {

    @synlist = ( '00', '06', '12', '18') ;


   } else {
    @synlist = ($syntime) ;

    }

     print "synlist = $synlist\n";



 $synlist_len=@synlist;
if ($synlist_len > 0){ 
$i=0;
   while  ($i < $synlist_len ){

     $synhour = $synlist[$i] ;

    print "synhour = $synhour\n";


#   gmi_L1CR.20140918.t12z.bufr
            $daily_bufr=token_resolve("$template_BUFR",$process_date,$synhour );
    print " POSCHITALI daily_bufr = $daily_bufr\n";


    if ( $synhour eq '00') {
          $knum = $k00;
          $clist =$list00;
          @lst_file = @lst_file00 ;
    }

    if ( $synhour eq '06') {
          $knum = $k06;
          $clist =$list06;
          @lst_file = @lst_file06 ;
    }


    if ( $synhour eq '12') {
          $knum = $k12;
          $clist =$list12;
          @lst_file = @lst_file12 ;
    }


    if ( $synhour eq '18') {
          $knum = $k18;
          $clist =$list18;
          @lst_file = @lst_file18 ;
    }

#            print "NADO   knum=    is $knum\n";
        $nol = 0;


          if ( $knum > $nol ) {
$inum=0;

         while  ($inum < $knum ){

#            $need_file =  $lst_file[$inum] ;
#        $rc=system("gmi1cr_bufr.x -d $process_date -t $synhour -f  $need_file");

         $rc=system("gmi1cr_bufr.x -d $process_date -t $synhour -f  $lst_file[$inum]");

 


if ($rc != 0) {

      print " WARNING:error running gmi1cr_bufr.x  for  $process_date synoptic $synhour. \n";
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: WARNING :error running gmi1cr_bufr.x for  $process_date synoptic $synhour"});
          $archive_err ++;
         recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);


              if ( ! $opt_f ) {
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: error running gmi1cr_bufr.x  for  $process_date synoptic $synhour"});

     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error: error running gmi1cr_bufr.x data for date: $process_date and for time: $synhour .  This error occurred while processing GMI  data .";
                              }

             }




         $inum = $inum + 1; 
#   end  for number of files 
          }
#  if  knum >0
      } 
       else {
    
      print " WARNING: NO input GMI data for  $process_date synoptic $synhour. \n";
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: WARNING :NO input GMI data for  $process_date synoptic $synhour"});
          $archive_err ++;
         recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);

              if ( ! $opt_f ) {
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: NO input GMI data for  $process_date synoptic $synhour"});

     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error: NO input GMI data for date: $process_date and for time: $synhour .  This error occurred while processing GMI  data .";
              }

       }
   

#----------------------------------------------------------------------------- 
# BEGIN  copy and archive daily_bufr    copy and archive daily_bufr  copy and archive daily_bufr
#   for $process_date , $synhour
#-----------------------------------------------------------------------------

# From Jianjun
#  if ( $size == 1654 || $size == 1697 || $size == 1840 || $size == 1883 || $size == 1926 ) 

       $min_size = 2000 ;
#    if (! -e "$daily_bufr" || -z "$daily_bufr") 
              $filesize = -s $daily_bufr ;
     if (! -e "$daily_bufr" || -z "$daily_bufr" ||  "$filesize" < $min_size ) {
        print "$daily_bufr DOES NOT EXIST\n";

        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: WARNING: $daily_bufr was NOT created,or is zero bytes,or size less than min size.  There may be no data available for date: $process_date and for times: $synhour .  This error occurred while processing GMI data ."});

     recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);

              if ( ! $opt_f ) {    
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: Error: $daily_bufr was NOT created,or is zero bytes,or size less than min size.  There may be no data available for date: $process_date and for times: $synhour .  This error occurred while processing GMI data ."});

     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
        die "Error: $daily_bufr was NOT created or is small size while running for GMI . ";
              }

    }
    else {
      print "$daily_bufr exists.  Will copy\n";
#  Copy into STAGE_DIR
    $rc=system("cp $daily_bufr $GMI_STAGE_DIR");

      if ($rc != 0) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: error copying daily_bufr file $daily_bufr."});
     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
        die "error copying daily_bufr file $daily_bufr.";
      }

# ---------------------------------------------
# Archive    $daily_bufr
#---------------------------------------------------


#     $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/gmi/bufr/Y%y4/M%m2/",$process_date);
#     $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'bufr', "$process_date",
#         "$archive_dir", "$daily_bufr", { 'verbose' => "$verbose" ,'exp_path' => "1" } );


      $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'bufr', "$process_date",
          "${GMI_ARCHIVE_LOC}", "$daily_bufr", { 'verbose' => "$verbose"  } );



      if ($rc != 1) {
      print " WARNING: Could not archive $daily_bufr .  \n";
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: WARNING :could not archive $daily_bufr"});
          $archive_err ++;
         recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);

              if ( ! $opt_f ) {
          err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: could not archive $daily_bufr"});

     recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
          die "Error: could not archive $daily_bufr  for date: $process_date and for time: $synhour .  This error occurred while processing GMI  data .";
           }
      }

#-----------------------------------------------------------------------------
# END  copy and archive daily_bufr    copy and archive daily_bufr  copy and archive daily_bufr
#-----------------------------------------------------------------------------

#   end  for else daily_bufr exist
  }



###############################

#   end  for foreach $synhour
          $i=$i+1;
    }
#  if synlist > 0
 }

########################
# create TAR files from raw GMI 1B nd 1CR files for $process date and 
# previous (if needed) date  and archive them
########################

#    Get  filename_beg_1CR and filename_beg_1B  first
   $GMI_ACQUIRE_PATH=$template_path_1B;
@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

 ${filename_beg_1B} = substr( ${GMI_ACQUIRE_FILE},0,12) ;

$template_path_1CR=$GMI_ACQUIRE_PATH_1CR ;

   $GMI_ACQUIRE_PATH=$template_path_1CR;
@list = split('/',${GMI_ACQUIRE_PATH});
$list_len = @list;
${GMI_ACQUIRE_FILE} = $list[$list_len-1];

 ${filename_beg_1CR} = substr( ${GMI_ACQUIRE_FILE},0,16) ;


          $nol = 0;

########################
# create TAR files from raw GMI 1B  files  for $process date and previous(if needed)
# date  and archive them
########################
#  01234567890123456789012345
#  1B.GPM.GMI.TB2014.20140416-S113146-E113644.V01D.RT-H5
#  1C-R.GPM.GMI.XCAL2014-N.20140416-S113146-E113644.V01D.RT-H5

# 03/17/2016  New name for HDF files
#  1C-R.GPM.GMI.XCAL2015-C.20160308-S235640-E000138.V04A.RT-H5   NEW HDF NAME
#  1B.GPM.GMI.TB2015.20160313-S000140-E000638.V04A.RT-H5   - NEW HDF Name  



#     $tar_name_1B=substr( ${filename_today_1B},0,18) ;
#     $tar_name_1CR=substr( ${filename_today_1CR},0,24) ;
#---------------------------------------------------------
#-----------------------------------------------------------
#  $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/gmi/native/1B/Y%y4/M%m2/",$process_date);


   get_native (${filename_beg_1B},$process_date) ;
      $len1B_today = $real_len_full ;
          print " len1B_today =  $real_len_full\n";

       if ($len1B_today >0) {
        ${need_name_1B} = $full_list[$nol];
 $tar_name_1B=substr( ${need_name_1B},0,18) ;


   $tar_name  =token_resolve("${tar_name_1B}%y4%m2%d2.he5.tar",$process_date);

                 print " iCHTORAT archive_dir = $archive_dir     \n";
                print " iCHTORAT tar_name = $tar_name     \n";
                print " iCHTORAT filename_today_1B = ${filename_today_1B}     \n";

#  TAR  Input Raw files for current day
   $rc=system(" tar -cvf ${tar_name} ${filename_today_1B} ");
#  Archive TAr file for current day
   $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'native', "$process_date",
          "$GMI_ARCHIVE_LOC", "$tar_name", { 'subtype' => "1B", 'verbose' => "$verbose"  } );

    if ($rc != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name  while running for GMI "});
        print "WARNING: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);
            if ( ! $opt_f ) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name while running for GMI"});
        print "ERROR: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);


          die "Error: could not archive $tar_name  for date: $process_date .  This error occurred while processing GMI data .";
           }
    }

# End  if ($len1B_today >0)
  }


########################
# create TAR files from raw  GMI 1CR  files for $process date and previous  date 
#  and archive them
########################

#  $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/gmi/native/1CR/Y%y4/M%m2/",$process_date);
# 1C-R.GPM.GMI.XCAL2014-N.20140924-S085646-E090144.V03B.RT-H5 
#  $tar_name     =token_resolve("1C-R.GPM.GMI.XCAL_%y4m%m2%d2.he5.tar",$process_date);

 get_native (${filename_beg_1CR},$process_date) ;
      $len1CR_today = $real_len_full ;
         print " len1CR_today =  $real_len_full\n";
    if ($len1CR_today >0) {
        ${need_name_1CR} = $full_list[$nol];
 $tar_name_1CR=substr( ${need_name_1CR},0,24) ;


   $tar_name =token_resolve("${tar_name_1CR}%y4%m2%d2.he5.tar",$process_date);

                 print " iCHTORAT archive_dir = $archive_dir     \n";
                print " iCHTORAT tar_name = $tar_name     \n";
                print " iCHTORAT filename_today_1CR = ${filename_today_1CR}     \n";

#  TAR  Input Raw files for current day
   $rc=system(" tar -cvf ${tar_name} ${filename_today_1CR} ");
#  Archive TAr file for current day
   $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'native', "$process_date",
          "$GMI_ARCHIVE_LOC", "$tar_name", {  'subtype' => "1CR",'verbose' => "$verbose"  } );

    if ($rc != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name  while running for GMI "});
        print "WARNING: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);
            if ( ! $opt_f ) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name while running for GMI"});
        print "ERROR: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);


          die "Error: could not archive $tar_name  for date: $process_date .  This error occurred while processing GMI data .";
           }
    }

#      if ($len1CR_today >0)
  }

#  444444444444444444444444444444444444444444444444444444444444444
     if ( ($syntime eq '00') or ($syntime eq 'hd'))  {
#  If we have to get 0z , or the whole day ( syntime = '00', or syntime = 'hd')
#  then we have files from the  previous day to create  TAR file for yesterday 
#    $tar_name_1B=substr( ${filename_yesterday_1B},0,18) ;
#     $tar_name_1CR=substr( ${filename_yesterday_1CR},0,24) ;


##################################
#  Tar files fot 1B for yesterday
##################################

#   $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/gmi/native/1B/Y%y4/M%m2/",$process_date_m1);
#  $tar_name     =token_resolve("1B.GMI-Aura_L2-OMTO3_%y4m%m2%d2.he5.tar",$process_date_m1);


 get_native (${filename_beg_1B},$process_date_m1) ;
      $len1B_yesterday = $real_len_full ;
         print " len1B_yesterday =  $real_len_full\n";
     if ($len1B_yesterday >0) {
        ${need_name_1B} = $full_list[$nol];
 $tar_name_1B=substr( ${need_name_1B},0,18) ;


 $tar_name =token_resolve("${tar_name_1B}%y4%m2%d2.he5.tar",$process_date_m1);

                 print " aCHTORAT archive_dir = $archive_dir     \n";
                print " aCHTORAT tar_name = $tar_name     \n";
                print " aCHTORAT filename_yesterday_1B = ${filename_yesterday_1B}     \n";


#  TAR  Input Raw files for previous day
   $rc=system(" tar -cvf ${tar_name} ${filename_yesterday_1B} ");
#  Archive  TAR  file for previous day
   $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'native', "$process_date_m1",
          "$GMI_ARCHIVE_LOC", "$tar_name", { 'subtype' => "1B", 'verbose' => "$verbose"  } );

   if ($rc != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name  while running for GMI "});
        print "WARNING: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);
            if ( ! $opt_f ) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name while running for GMI"});
        print "ERROR: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);


          die "Error: could not archive $tar_name  for date: $process_date_m1 .  This error occurred while processing GMI data .";
           }
    }
#   if ($len1B_yesterday >0)
  }


#####################################
#  Tar files  1CR for yesterday
#####################################
#   $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/gmi/native/1CR/Y%y4/M%m2/",$process_date_m1);

 get_native (${filename_beg_1CR},$process_date_m1) ;
      $len1CR_yesterday = $real_len_full ;
 print " len1CR_yesterday =  $real_len_full\n";
       if ($len1CR_yesterday >0) {
        ${need_name_1CR} = $full_list[$nol];
 $tar_name_1CR=substr( ${need_name_1CR},0,24) ;



 $tar_name=token_resolve("${tar_name_1CR}%y4%m2%d2.he5.tar",$process_date_m1);

                 print " aCHTORAT archive_dir = $archive_dir     \n";
                print " aCHTORAT tar_name = $tar_name     \n";
                print " aCHTORAT filename_yesterday_1CR = ${filename_yesterday_1CR}     \n";


#  TAR  Input Raw files for previous day
   $rc=system(" tar -cvf ${tar_name} ${filename_yesterday_1CR} ");
#  Archive  TAR  file for previous day
   $rc = gen_archive ( "$env", "$prep_ID", 'gmi', 'native', "$process_date_m1",
          "$GMI_ARCHIVE_LOC", "$tar_name", { 'subtype' => "1CR", 'verbose' => "$verbose"  } );

   if ($rc != 1) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name  while running for GMI "});
        print "WARNING: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, WARNING, $tab_argv, $sched_dir, $sched_sts_fl);
            if ( ! $opt_f ) {
        err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
                 {'err_desc' => "${err_pref}: could not archive $tar_name while running for GMI"});
        print "ERROR: could not archive $tar_name \n";
        $archive_err ++;
          recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);


          die "Error: could not archive $tar_name  for date: $process_date_m1 .  This error occurred while processing GMI data .";
             }

      }

# END       if ($len1CR_yesterday >0)
  }

# END     if ( ($syntime eq '00') or ($syntime eq 'hd'))
 }



########################
# Rename output listings
########################

if ( $opt_O ) {

    print "THERE ARE  OPTION O.  Will copy Listing\n";

    unlink<"$listing_file_gz">;
    system ( "gzip -c $listing_file > $listing_file_gz" );

    $archive_dir = token_resolve("${GMI_ARCHIVE_LOC}/listing/Y%y4/M%m2/",$process_date);
    $rc=gen_archive ( "$env","$prep_ID",'gmi','listing', "$process_date",
                     "$GMI_ARCHIVE_LOC", $listing_file_gz,
                     { 'remote_name' => "gmi_${prep_ID}.$err_time.listing.gz",
                       'delete'      => "1",
                       'verbose'     => "1" } );

    if ( $rc != 1 ) {
	err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1", {'err_desc' => "${err_pref}: could not archive listing file gmi_${prep_ID}.$err_time.listing.gz"});
	print "WARNING: could not archive listing file gmi_${prep_ID}.$err_time.listing.gz \n";
	$archive_err ++;
    }
    
} 
if ( $archive_err ) {
    err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
	     {'err_desc' => "${err_pref}: get_gmi.pl: exiting with errors"});
}else{
    err_log (0, "get_gmi.pl", "$err_time","$prep_ID","-1",
	     {'err_desc' => "${err_pref}: get_gmi.pl: exiting normally"});
}
if ( $opt_O ) {
    system ("mv $listing_file $opt_O/gmi_${prep_ID}.${err_time}.listing");
}
############################
# Clean up working directory 
############################

  $rc=system("/bin/rm -rf ${gmi_work}");


   if ($rc != 0) {
    err_log (4, "get_gmi.pl", "$err_time","$prep_ID","-1",
             {'err_desc' => "${err_pref}: WARNING: could not remove ${gmi_work}"});
    print "WARNING: could not remove ${gmi_work}\n";
    }




   recd_state( $fl_name, "COMPLETE", $tab_argv, $sched_dir, $sched_sts_fl );
exit 0;


