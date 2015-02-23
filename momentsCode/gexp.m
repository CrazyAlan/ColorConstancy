function imout = gexp(imin)

a= 10;

% GLOG:  imout = a*( imin.^(1/a) - 1 ) ;

% So, inverse is:
imout = ( imin/a+1.0 ).^a;
