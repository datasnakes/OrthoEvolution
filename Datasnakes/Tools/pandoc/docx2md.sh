#!/bin/bash

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
