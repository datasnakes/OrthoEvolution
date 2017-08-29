#!/bin/bash

# This script was created to convert a directory full
# of Word Document (.docx) files into .md equivalents. It uses
# pandoc to do the conversion.
#
# 1. Install pandoc from http://johnmacfarlane.net/pandoc/
# 2. Copy this script into the directory containing the .md files
# 3. Ensure that the script has execute permissions
# 4. Run the script

# By default this will keep the original .md file
# This script also removes white space from file names.

# White space removal
FILES1=*.docx
for f in $FILES1
do
  mv "$f" "${f//[[:space:]]}"
  echo "Removed whitespace from $f"
done

# DOCX to MD conversion
FILES2=*.docx
for f in $FILES2
do
  # extension="${f##*.}"
  filename="${f%.*}"
  echo "Converting $f to $filename.md"
  `pandoc "${f}" -t markdown -o $filename.md`
  # uncomment this line to delete the source file.
  # rm $f
done

# Recreates the original file
for f in $FILES1
do
