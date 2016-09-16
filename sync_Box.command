#!/usr/bin/env bash

folder="/Users/lucasnell/Box Sync/ZooEnt_540_2016"

# Change directory to this file's location
cd `dirname $0`

rsync -razP --delete --exclude "Readings" "$folder" ./box