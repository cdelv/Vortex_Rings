#!/bin/sh

darcsroot()
{
    d=`pwd`
    while ! test -d _darcs; do
	cd ..
    done
    pwd
    cd "$d"
}

showfiles()
{
    if ! darcs show files > /dev/null 2>&1; then
	ls *.$1 2> /dev/null | sed 's/\(.*\)/	\1 \\/g'
    else
	ROOT=`darcsroot`
	DIR=`echo $PWD | sed "s|$ROOT|.|"`
	darcs show files | grep '\'$DIR'/[^/ ]*\.'$1'$' | \
	    sed -e 's|'$DIR'/\(.*\)|	\1 \\|g'
    fi
}

showpages()
{
    if ! darcs show files > /dev/null 2>&1; then
	files=`ls *$1 2> /dev/null`
    else
	ROOT=`darcsroot`
	DIR=`echo $PWD | sed "s|$ROOT|.|"`
	files=`darcs show files | grep '\'$DIR'/[^/ ]*'$1'$' | sed "s|$DIR/||g"`
    fi
    for f in $files; do
	case "$f" in
	    *.page) echo "$f" ;;
	    *.[chm]) $BASILISK/darcsit/pagemagic "$f" && echo "$f.page" ;;
	esac
    done | sed 's/\(.*\)/	\1 \\/g'
}

# join lines delimited by \\n characters
singleline()
{
    sed -e ':a
N
$!ba
s/\\\n//g
'
}

echo "updating Makefile.tests"
(
    echo "# Automatically generated using 'make Makefile.tests'"
    echo "# DO NOT EDIT, edit 'Makefile' instead"
    echo "ALLTESTS = \\"
    (showfiles c
     grep '^[a-zA-Z_0-9-]*\.*tst[ ]*:' Makefile | \
	 sed -n 's/\(^[a-zA-Z_0-9-]*\)\.*tst[ ]*:.*/	\1.c \\/p'
    ) | sort | uniq
    echo ""
    echo "plots: \\"
    showfiles plot | sed 's/\(.*\)\.plot/\1\/plot.png/g'
    showpages '.c' | sed 's/\(.*\)\.c/\1\/plots/g'
    echo ""
    echo "TESTS = \\"
    singleline < Makefile | 		                   \
	grep '^check:' | tr ' ' '\n' |                     \
	sed -n 's/[ 	]*\([a-zA-Z_0-9-]*\..tst\)[ 	]*/\1/p' |  \
	sed 's/\(.*\)\..tst/	\1.c \\/g'
    singleline < Makefile | 		                   \
	grep -v '^check:' | grep '^[^.]*:.*' | tr ' ' '\n' |     \
	sed -n 's/[ 	]*\([a-zA-Z_0-9-]*\.tst\)[ 	]*/\1/p' |    \
	sed 's/\(.*\)\.tst/	\1.c \\/g'
    echo ""
    echo "SPECIAL_TESTS = \\"
    sed -n 's/.*:[ 	]*\([a-zA-Z_0-9-]*\.ctst\)/\1/p' Makefile | \
	sort | uniq | sed 's/\(.*\)/	\1 \\/g'
    echo ""
    echo "ALLPAGES = \\"
    showpages ''
    echo ""
) > Makefile.tests

rm -f Makefile.deps
