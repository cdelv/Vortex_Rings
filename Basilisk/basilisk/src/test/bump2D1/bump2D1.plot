! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 400,400

set size ratio -1
unset key
unset xtics
unset ytics
unset border
unset colorbox
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

set multiplot layout 3,3 scale 1.75,1.75
splot 'eta-0'
splot 'eta-1'
splot 'eta-2'
splot 'eta-3'
splot 'eta-4'
splot 'eta-5'
splot 'eta-6'
splot 'eta-7'
splot 'eta-8'
unset multiplot

set output 'level.png'

set multiplot layout 3,3 scale 1.75,1.75
splot 'level-0'
splot 'level-1'
splot 'level-2'
splot 'level-3'
splot 'level-4'
splot 'level-5'
splot 'level-6'
splot 'level-7'
splot 'level-8'
unset multiplot

exit 0

set output 'pid.png'

set multiplot layout 3,3 scale 1.6,1.6
splot 'pid-0'
splot 'pid-1'
splot 'pid-2'
splot 'pid-3'
splot 'pid-4'
splot 'pid-5'
splot 'pid-6'
splot 'pid-7'
splot 'pid-8'
unset multiplot

! rm -f eta-? level-? pid-?
