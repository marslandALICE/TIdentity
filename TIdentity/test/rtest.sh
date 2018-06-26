#!/bin/bash

dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/DataTree.root
lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/LineShapes.root
sign=$1

./testIdenTest $dataTree $lineShapes  $sign
