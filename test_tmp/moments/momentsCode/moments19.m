function rowvec = moments19(im)
    R = mean(im(:,1));
    G = mean(im(:,2));
    B = mean(im(:,3));
    R2 = sqrt(mean(im(:,1).^2));
    G2 = sqrt(mean(im(:,2).^2));
    B2 = sqrt(mean(im(:,3).^2));
    RG = sqrt(mean(im(:,1).*im(:,2)));
    RB = sqrt(mean(im(:,1).*im(:,3)));
    GB = sqrt(mean(im(:,2).*im(:,3)));
    
    R3 = (mean(im(:,1).^3)).^1/3;
    G3 = (mean(im(:,2).^3)).^1/3;
    B3 = (mean(im(:,3).^3)).^1/3;
    R2G = (mean((im(:,1).^2).*im(:,2))).^1/3;
    R2B = (mean((im(:,1).^2).*im(:,3))).^1/3;
    G2R = (mean((im(:,2).^2).*im(:,1))).^1/3;
    G2B = (mean((im(:,2).^2).*im(:,3))).^1/3;
    B2R = (mean((im(:,3).^2).*im(:,1))).^1/3;
    B2G = (mean((im(:,3).^2).*im(:,2))).^1/3;
    RGB = (mean(im(:,1).*im(:,2).*im(:,3))).^1/3;
    
    
    
    rowvec = [R,G,B,R2,G2,B2,RG,RB,GB,R3,G3,B3,R2G,R2B,G2R,G2B,B2R,B2G,RGB];
return