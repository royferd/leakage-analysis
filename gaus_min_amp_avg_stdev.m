function gaus = gaus_min_amp_avg_stdev(x,xdata,ydata)

A = x(1);

%B = xavg;
B = x(2);


%C = xstd;
C = x(3);

gaus = sum( ( ydata - A^2*exp( -( ( xdata - B )/(sqrt(2)*C) ).^2 ) ).^2 );