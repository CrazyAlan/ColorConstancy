function imout = glog(imin)

a= 10;

imout=imin;

temp = (imin>=0);
imout(temp) = ( imin(temp).^(1/a) - 1 ) ./ (1/a);
imout(~temp) = - ( (-imin(~temp)).^(1/a) - 1 ) ./ (1/a);


