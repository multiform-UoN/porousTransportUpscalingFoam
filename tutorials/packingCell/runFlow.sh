#!/bin/bash



. $WM_PROJECT_DIR/bin/tools/RunFunctions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions


runApplication simpleFoam
foamEndTime=$(foamListTimes | tail -1)
if [ -d $foamEndTime ]; then
   mkdir 0.end
   mv $foamEndTime/U 0/.
   rm -r $foamEndTime
fi
