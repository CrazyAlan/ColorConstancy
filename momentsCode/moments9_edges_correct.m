function rowvec = moments9_edges(dx,dy)
    Rx= mean(dx(:,1)); % scalar
    Gx= mean(dx(:,2));
    Bx= mean(dx(:,3));
    R2x= sqrt(mean(dx(:,1).^2));
    G2x= sqrt(mean(dx(:,2).^2));
    B2x= sqrt(mean(dx(:,3).^2));
%     RGx= sqrt(mean( abs(dx(:,1).*dx(:,2))) );
%     RBx= sqrt(mean( abs(dx(:,1).*dx(:,3))) );
%     GBx= sqrt(mean( abs(dx(:,2).*dx(:,3))) );
    tt=mean( (dx(:,1).*dx(:,2)) );  RGx= sign(tt)*sqrt(abs(tt));
    tt=mean( (dx(:,1).*dx(:,3)) );  RBx= sign(tt)*sqrt(abs(tt));
    tt=mean( (dx(:,2).*dx(:,3)) );  GBx= sign(tt)*sqrt(abs(tt));
    rowvecx = [Rx,Gx,Bx,R2x,G2x,B2x,RGx,RBx,GBx];
    %
    Ry= mean(dy(:,1));
    Gy= mean(dy(:,2));
    By= mean(dy(:,3));
    R2y= sqrt(mean(dy(:,1).^2));
    G2y= sqrt(mean(dy(:,2).^2));
    B2y= sqrt(mean(dy(:,3).^2));
%     RGy= sqrt(mean( abs(dy(:,1).*dy(:,2))) );
%     RBy= sqrt(mean( abs(dy(:,1).*dy(:,3))) );
%     GBy= sqrt(mean( abs(dy(:,2).*dy(:,3))) );
    tt=mean( (dy(:,1).*dy(:,2)) );  RGy= sign(tt)*sqrt(abs(tt));
    tt=mean( (dy(:,1).*dy(:,3)) );  RBy= sign(tt)*sqrt(abs(tt));
    tt=mean( (dy(:,2).*dy(:,3)) );  GBy= sign(tt)*sqrt(abs(tt));
    rowvecy = [Ry,Gy,By,R2y,G2y,B2y,RGy,RBy,GBy];
%
rowvec=[rowvecx,rowvecy]; % 1 x 18
return