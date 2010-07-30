pause 0.002
t=t+dt
#set view 80, (t/dt + 15)%360
replot;
if (t+dt<=T) load 'timeStepping.gnu'

