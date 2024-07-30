###############################################################
##      Script to Transform Reduze Output to FORM Input      ##
##                     01May15, IMSc                         ##
##                        Taushif                            ##
###############################################################

#!/bin/bash

echo 'Script to Transform Reduze Output to FORM Input for 2 loop'
echo 'Enter the name of the .inc file (2-loop)!!'
read input
echo 'Enter the desired output name !!'
read output

#-------------------------------------
# Removing the bracket around DiaMatch 
#-------------------------------------

sed -r 's|DiaMatch\(([0-9]+)\)|DiaMatch\1|g' $input > 1.d

#-----------
# Removing *
#-----------

sed -r 's|\*||g' 1.d > 2.d

#-----------
# Removing =
#-----------

sed -r 's|\=||g' 2.d > 3.d

#-------------
# Removing [],
#-------------

sed -r 's|\[\]\, ||g' 3.d > 4.d

#---------------------------------------------
# Removing , [] The last one appears like this
#---------------------------------------------

sed -r 's|\, \[\]||g' 4.d > 5.d

#------------------------------------------
# Converting Sector(A1, x, y) to multiplyA1
#------------------------------------------

sed -r 's|Sector\(([^ ]*)\, ([^ ]*) ([^ ]*)\)|multiply\1\;|g' 5.d > 6.d

#-------------
# Removing id 
#-------------

sed -r 's|id ||g' 6.d > 7.d

#-------------------------------------
# Replacing double space by single one
#-------------------------------------

sed -r 's|  | |g' 7.d > 8.d

#----------------------------------------------------
# Copying DiaMatch'i' to MatchDia'i' in the next line
#----------------------------------------------------

sed -r 's|DiaMatch([0-9]+) ([^ ]*) ([^ ]*) ([^ ]*)\;|DiaMatch\1 \2 \3 \4 \;\nMatchDia\1\n|g' 8.d > 9.d

#------------------------------------------------------------------------
# Keep the DiaMatch'i' in the same line and put the rest in the next line
#------------------------------------------------------------------------

sed -r 's|DiaMatch([0-9]+) ([^ ]*) ([^ ]*) ([^ ]*)|DiaMatch\1 \n \2 \3 \4|g' 9.d > 10.d

#-------------------------------------------------
# Change the DiaMatch'i' to #procedure DiaMatch'i'
#-------------------------------------------------

sed -r 's|DiaMatch([0-9]+)|\#procedure DiaMatch\1|g' 10.d > 11.d

#----------------------------------------------------
# Change the MatchDia'i' to #endprocedure DiaMatch'i'
#----------------------------------------------------

sed -r 's|MatchDia([0-9]+)|\#endprocedure DiaMatch\1|g' 11.d > 12.d

#-------------------------------------------------
# replace multiply by multiply followed by a space
#-------------------------------------------------

sed -r 's|multiply([^ ]*)|multiply \1|g' 12.d > 13.d

#-----------------------------------
# Replace Shift by multiply replace_
#-----------------------------------

sed -r 's|Shift|multiply replace\_|g' 13.d > $output

#------------------
# Removing .d files 
#------------------

rm *.d
echo 'Enjoy!!'
