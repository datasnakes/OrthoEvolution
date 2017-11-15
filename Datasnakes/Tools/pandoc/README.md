# Pandoc Documentation
Use the `docx2md.sh` script to convert .docx files to .md (markdown) format.
The shell script uses pandoc to convert the files.

## Dependencies
[Pandoc](http://johnmacfarlane.net/pandoc/) must be installed.

## Setup
**On Linux/Debian**

Make the script executable. Then run it.
1. `chmod +x docx2md.sh`
2. `./docx2md.sh`

## Examples
In addition to a .sh/bash script to use with Pandoc, we've used [pypandoc]() to
create a class that allows the conversion of documents.

### Convert docx to markdown
```python
```