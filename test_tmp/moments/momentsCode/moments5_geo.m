function rowvec = moments5_geo(im)
    R = mygeo_glog(im(:,1));
    G = mygeo_glog(im(:,2));
%     B = mygeo_glog(im(:,2));
    R2 = sqrt(mygeo_glog(im(:,1).^2));
    G2 = sqrt(mygeo_glog(im(:,2).^2));
%     B2 = sqrt(mygeo_glog(im(:,3).^2));
    RG = sqrt(mygeo_glog(im(:,1).*im(:,2)));
%     RB = sqrt(mygeo_glog(im(:,1).*im(:,3)));
%     GB = sqrt(mygeo_glog(im(:,2).*im(:,3)));
%     rowvec = [R,G,B,R2,G2,B2,RG,RB,GB];
    rowvec = [R,G,R2,G2,RG];
return