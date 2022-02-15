#!/bin/bash
date

if [ $# -eq 1 ]

type=$1 #Inner, Outer

then
  echo ${type}
  rm ${type}Sensors.list

  for FILE in `cat FstModules.list`
  do
       echo "${FILE}"
       ls /Users/gavinwilks/Workspace/ForwardSiliconTracker/FstSurvey/FstSurvey/SensorSurveyMatrices/data/${type}Sensors/*${FILE}* >> ${type}Sensors.list
  done
fi

