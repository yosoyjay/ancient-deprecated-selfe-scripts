#!/bin/bash
if [ ! -z "$STY" ] ; then
  echo Not a good idea to run this from inside screen
    exit 0
fi
PTS=`ps -ef | grep screen | grep -v grep | grep $LOGNAME | awk '{print $6}'`
if [ ! -z "$PTS" ] ; then
  PS=`ps -ef | grep $PTS | grep bash | awk '{print $2}'`
  if [ ! -z "$PS" ] ; then
  	echo killing $PS
    kill -9 $PS
    screen -list
  fi
fi
