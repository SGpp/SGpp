#set hidden3d

dt=16384; ## step through lines by this much
T=100*dt; ## go until this line reached
t=99*dt;   ## start at line zero

set xrange [0:1]
set yrange [0:1]  

scale = 0.2
plot "out.dat" every ::t::t+dt using 1:2:($3*scale):($4*scale) with vector title ""

#load "timeStepping.gnu"
