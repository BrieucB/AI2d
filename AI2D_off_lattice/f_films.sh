#!/bin/bash

make clean ; make

cd films_nei_area_unit_v1

for D in 0.05 0.1 0.15 0.2
do
    cd D$D
    printf -v D0 "%0.2f" $D
    cd run_rho0.95_v1.77_D${D0}_beta2.00/snaps_prof
    ffmpeg -i snap_%05d.png -c:v libx264 -pix_fmt yuv420p ../../../rho0.95_v1.77_D${D}_beta2.00.mp4
    cd ../../..
done
