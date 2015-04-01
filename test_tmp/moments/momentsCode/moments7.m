function rowvec = moments9(im)
    R = mean(im(:,1));
    G = mean(im(:,2));
    B = mean(im(:,3));
    R2 = sqrt(mean(im(:,1).^2));
    G2 = sqrt(mean(im(:,2).^2));
    B2 = sqrt(mean(im(:,3).^2));
    RG = sqrt(mean(im(:,1).*im(:,2)));
    RB = sqrt(mean(im(:,1).*im(:,3)));
    GB = sqrt(mean(im(:,2).*im(:,3)));
    rowvec = [R,G,B,R2,G2,B2,RG,RB,GB];
return