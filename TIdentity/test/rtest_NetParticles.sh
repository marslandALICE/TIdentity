#!/bin/bash

dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/DataTree_Net.root
lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/LineShapes_Net.root
sign=$1

./testIden_NetParticles $dataTree $lineShapes  $sign
