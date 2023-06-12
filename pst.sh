#!/bin/bash

echo "Post-processing script started..."
python3 pltmsh.py
ffmpeg -i pics/v_jpg/v%d.jpg -vcodec mpeg4 v.avi
echo "Post-processing script complete..."
