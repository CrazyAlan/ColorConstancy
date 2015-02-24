function zeta = zetaIm(roke,im)
    imchrom = makechrom3vec(im); % n by 3
    phik = imchrom./(repmat(roke,size(imchrom,1),1));
    phik = glog(phik);
    zeta = -phik*roke';
return