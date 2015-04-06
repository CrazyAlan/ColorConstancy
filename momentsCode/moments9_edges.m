function rowvec = moments9_edges(dx,dy)
    Rx= mean(dx(:,1)); % scalar
    Gx= mean(dx(:,2));
    Bx= mean(dx(:,3));
    R2x= sqrt(mean(dx(:,1).^2));
    G2x= sqrt(mean(dx(:,2).^2));
    B2x= sqrt(mean(dx(:,3).^2));
    RGx= sqrt(abs(mean(dx(:,1).*dx(:,2))));
    RBx= sqrt(abs(mean(dx(:,1).*dx(:,3))));
    GBx= sqrt(abs(mean(dx(:,2).*dx(:,3))));
    rowvecx = [Rx,Gx,Bx,R2x,G2x,B2x,RGx,RBx,GBx];
    %
    Ry= mean(dy(:,1));
    Gy= mean(dy(:,2));
    By= mean(dy(:,3));
    R2y= sqrt(mean(dy(:,1).^2));
    G2y= sqrt(mean(dy(:,2).^2));
    B2y= sqrt(mean(dy(:,3).^2));
    RGy= sqrt(abs(mean(dy(:,1).*dy(:,2))));
    RBy= sqrt(abs(mean(dy(:,1).*dy(:,3))));
    GBy= sqrt(abs(mean(dy(:,2).*dy(:,3))));
    rowvecy = [Ry,Gy,By,R2y,G2y,B2y,RGy,RBy,GBy];
    %
    rowvec=[rowvecx,rowvecy]; % 1 x 18
return
