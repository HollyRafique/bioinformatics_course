#!/usr/bin/env bash
#empty dirs in repo
mkdir -p analysis docs data

#adding READMEs to each dir
for my_directory in scripts analysis docs data; do
	touch ${my_directory}/README.md
	echo "# ${my_directory}" >> ${my_directory}/README.md
done

