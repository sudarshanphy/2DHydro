#!/bin/bash

#code to use the correct make file
DIR=$1
COMMAND=$2
cp ./sites/$DIR/Makefile .
make $COMMAND
