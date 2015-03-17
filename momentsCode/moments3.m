function rowvec = moments3(im)
    R = mean(im(:,1));
    G = mean(im(:,2));
    B = mean(im(:,3));
    
    rowvec = [R,G,B];
return