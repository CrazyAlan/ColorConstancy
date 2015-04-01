function chrom=makechrom3vec(mat)
% mat is Nx3
s = size(mat);
len = s(1);
mat = double(mat);
denom = mat(:,1)+mat(:,2)+mat(:,3);
back = denom==0;
denom(back)=999;
chrom = zeros(len,3);
for k=1:3
  chrom(:,k) = mat(:,k)./denom;
end
if (sum(back)~=0)
  chrom(back,:) = 0;
end
