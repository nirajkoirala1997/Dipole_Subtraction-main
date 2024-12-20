#!/bin/bash

filename="compare3.f"


home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)

cd summary/compare
gfortran $filename && ./a.out
rm a.out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          General Instructions for the compare.f Use                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#      -----------------------------------------------------------------------------------------------
#	compare1.f	is the older version of compare2.f can read data without error                 ‼️
#      -----------------------------------------------------------------------------------------------

#      -----------------------------------------------------------------------------------------------
#       compare2.f	has ability to display data from files directed by output_files.dat with error ✅
#      -----------------------------------------------------------------------------------------------

#      -----------------------------------------------------------------------------------------------
#	compare3.f 	can combine two files or by editing one can take ratio, difference and so on   ✅
#      -----------------------------------------------------------------------------------------------

#      -----------------------------------------------------------------------------------------------
#	compare4.f 	is the older version of compare5.f enter three files manually                  ‼️
#      -----------------------------------------------------------------------------------------------

#      -----------------------------------------------------------------------------------------------
#	compare5.f	this can combine plus regular and delta to userdefined output combined PK file ✅
#      -----------------------------------------------------------------------------------------------

#      -----------------------------------------------------------------------------------------------
#	compare6.f	this is the combination of compare2.f and compare5.f                           ✅
#      -----------------------------------------------------------------------------------------------


#      -----------------------------------------------------------------------------------------------
#	compare7.f  combinati  of compare2.f and compare5.f and can display PK indiv data + delta reg  ✅
#      -----------------------------------------------------------------------------------------------

#	✅ Working
#	‼️  No Longer updating
