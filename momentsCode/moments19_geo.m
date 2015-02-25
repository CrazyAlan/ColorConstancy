function rowvec = moments19_geo(im)
    R = mygeo_glog(im(:,1));
    G = mygeo_glog(im(:,2));
    B = mygeo_glog(im(:,3));
    R2 = sqrt(mygeo_glog(im(:,1).^2));
    G2 = sqrt(mygeo_glog(im(:,2).^2));
    B2 = sqrt(mygeo_glog(im(:,3).^2));
    RG = sqrt(mygeo_glog(im(:,1).*im(:,2)));
    RB = sqrt(mygeo_glog(im(:,1).*im(:,3)));
    GB = sqrt(mygeo_glog(im(:,2).*im(:,3)));
    
    R3 = (mygeo_glog(im(:,1).^3)).^1/3;
    G3 = (mygeo_glog(im(:,2).^3)).^1/3;
    B3 = (mygeo_glog(im(:,3).^3)).^1/3;
    R2G = (mygeo_glog((im(:,1).^2).*im(:,2))).^1/3;
    R2B = (mygeo_glog((im(:,1).^2).*im(:,3))).^1/3;
    G2R = (mygeo_glog((im(:,2).^2).*im(:,1))).^1/3;
    G2B = (mygeo_glog((im(:,2).^2).*im(:,3))).^1/3;
    B2R = (mygeo_glog((im(:,3).^2).*im(:,1))).^1/3;
    B2G = (mygeo_glog((im(:,3).^2).*im(:,2))).^1/3;
    RGB = (mygeo_glog(im(:,1).*im(:,2).im(:,3))).^1/3;
    
    
    
    rowvec = [R,G,B,R2,G2,B2,RG,RB,GB,R3,G3,B3,R2G,R2B,G2R,G2B,B2R,B2G,RGB];
return