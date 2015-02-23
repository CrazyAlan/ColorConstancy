function g = mygeo_glog(x)
% exp(1/sum(bright)*sum(log(imblurred(bright,:))))
%g = gexp( 1/length(x)*sum(glog(x)) );
mask=(x==0);
x=x(~mask);
g = gexp( 1/length(x)*sum(glog(x)) );