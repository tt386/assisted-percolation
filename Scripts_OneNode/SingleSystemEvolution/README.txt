The idea of this directory is to create a single repeat of a system, and plot the evolution of such a system as a range of images that can then be collated into a video using 

ffmpeg -r 10 -pattern_type glob -i 'AField_*.png' -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf scale=-2:1240 SystemEvolution.mp4


