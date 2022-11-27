#!/bin/bash
# This script needs one argument
# argc -> $#   ;    argv[1] = $1
if [ $# -ne 1 ]
then
    echo "Pass one argument"
    exit 0
fi

# add to git stage
git add .

# commit the changes
git commit -m "$1"

# push it to github server
git push
