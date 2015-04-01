function [rgbx,rgby] = makegradient(rgb)
%
[M,N,n3] = size(rgb);
if (M==1) || (N==1)
  disp('error in makegradient: must not be long-index image.');
  return;
end
tx=[-1,1];
ty=[-1;1];
%  freqz2([1,-1])
%
rgbx = zeros(M,N,n3);
rgby = zeros(M,N,n3);
for k=1:n3
 temp=conv2(rgb(:,:,k),fliplr(tx));
 rgbx(:,:,k)= temp(:,2:end);
 rgbx(:,end,k) = 0;
 temp=conv2(rgb(:,:,k),flipud(ty));
 rgby(:,:,k)= temp(2:end,:);
 rgby(end,:,k) = 0;
end
