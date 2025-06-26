#!/bin/bash


test=$1


t1="${test%.sorted.paf}"
t2="${t1##*to_}"

echo "${t2}"
