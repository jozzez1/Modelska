#!/bin/sh

CWD=`pwd`
DIR=`pwd`/$1

FPS=34

cd $DIR

mencoder -msgcolor "mf://*.jpg" -mf w=800:h=400:fps=${FPS}:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1.avi

cp $1.avi $CWD

cd $CWD
mplayer $1.avi

exit 0

