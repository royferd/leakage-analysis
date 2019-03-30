%function gaus = gaus_min_amp(x,xdata,ydata,xavg,xstd)
function gaus = gaus_min_amp(x,xdata,ydata,xavg,xstd)
%function gaus = gaus_min_amp(x,xdata,ydata)

A = x(1);

%B = -666.1133;
B = xavg;
%B = x(2);

%C = sqrt(2)*9.8417;
C = sqrt(2)*xstd;
%C = x(3);

%lambda = x(2);
%gaus = sum((ydata - A*exp(-lambda*tdata)).^2);
%gaus = sum((ydata - A*exp(-(xdata-xavg)^2/(2*xstd^2)))^2);
gaus = sum( ( ydata - A*exp( -( ( xdata - B )/C ).^2 ) ).^2 );