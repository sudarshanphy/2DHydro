#!/bin/bash

#if we need to loop over fields
fields="dens" 

# Loop over all fields
for field in $fields
do
ffmpeg -i dens_%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" dens.mp4
done

echo Finished loop over all fields.
