set hidden3d

t=1;   ## start at line zero
dt=10000; ## step through lines by this much
T=100*dt; ## go until this line reached


set zrange [0:9]  

splot "out.dat" every ::t::t+dt using 1:2:3 title "Time: ".t with lines

load "timeStepping.gnu"
