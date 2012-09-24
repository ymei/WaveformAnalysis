#!/bin/bash


# Convert init.scm into a C string constant so that it can be compiled into
# the executable

ofName=$1.h
echo "static const char init_scm_string[] = " > "$ofName"
sed -e 's/\"/\\"/g' \
    -e 's/^/"/g' \
    -e 's/$/\\n"/g' \
    $1 >> "$ofName"
echo ";" >> "$ofName"

