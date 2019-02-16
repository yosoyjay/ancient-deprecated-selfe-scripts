#!/bin/sh
#
# This script is used to call matlab functions from the command line.
#

matlab_exec=matlab
X="${1} ${2}"
echo ${X} > matlab_command_${2}.m
cat matlab_command${2}.m
${matlab_exec} -nodisplay -nosplash < matlab_command_${2}.m
rm matlab_command_${2}.m
