function [TIMESTEP,TEND] = wakepoints(n)


dAzi = 10;

RevTIMEENDAERO = 40;
RevSTREAMWISEWAKEPOINTS = 40;
RevFREESTRMWISEWAKEPOINTS = 40;

TIMEENDAERO = (RevTIMEENDAERO*360)/dAzi;
STREAMWISEWAKEPOINTS = (RevSTREAMWISEWAKEPOINTS*360)/dAzi;
FREESTRMWISEWAKEPOINTS = (RevFREESTRMWISEWAKEPOINTS*360)/dAzi;

TIMESTEP = (dAzi*pi/180)/(2*pi*n/60);
TEND = TIMEENDAERO*TIMESTEP;

end





