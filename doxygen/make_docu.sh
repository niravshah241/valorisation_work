#!/bin/bash

doxydir=$(find -wholename '*/Doxyfile.in')
doxydir=$(dirname $doxydir)
cd $doxydir
if [ ! -f ./Doxyfile.in ]; then
    echo "Error: make_docu.sh needs to be run in doxygen directory!"
    exit 1;
fi

if [ ! -f ./configuration ]; then
    echo "configuration file does not exist. Shall it be constructed? (y/n)"
    read answer
    [ $answer == "n" ] && exit 2;

    echo

    BASEDIR_GUESS=${PWD/\/doxygen//}
    echo "Please give me the base directory of the RBmatlab installation."
    echo "empy line for default = [ $BASEDIR_GUESS ]"
    read BASEDIR
    [ -z "$BASEDIR" ] && BASEDIR=$BASEDIR_GUESS

    echo 'define(`BASEDIR'\'', `'"$BASEDIR"\'')' > ./configuration
fi

m4 ./Doxyfile.in > ./Doxyfile

#m4 ./Doxyfile_latex.in > ./Doxyfile_latex

cd ..
doxygen doxygen/Doxyfile 2>&1 1>doxygen/doxygen.log | grep -v synupdate | grep -v docupdate | grep -v display | grep -v subsref | tee doxygen/doxygen.err
#doxygen doxygen/Doxyfile_latex 2> doxygen/doxygen_latex.err 1> doxygen/doxygen_latex.log

cd doxygen/html;

for i in *.html; do ../../bin/postprocess $i; done

#cd ../latex

#for i in *.tex; do ../../bin/postprocess $i; done
