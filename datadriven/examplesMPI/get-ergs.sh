#!/bin/bash
# grep for the results and remove the fluff around it - not exactly beautiful but it gets the job done
output="$(grep -r 'elapsed time:' $1 | sed 's/\b\posterrun-run.* \b//g' | sed 's/\//\ /g' | sed "s/\bBinary.*matches\b//g" | sed "s/\b$1\b//g" | sed '/^\s*$/d' | sed 's/.Nodes.//g' | sed 's/[s]//g')"
echo "$output"
# also output to file if one is specified
if ! [ -z "$2" ]
then
	echo "$output" > "$2"
fi
