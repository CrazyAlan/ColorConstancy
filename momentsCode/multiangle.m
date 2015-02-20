function angles = multiangle(actual,approx)
% actual and approx are both matrices.
% better be N x vec-dimension (say 3, or 7...)

ss1=size(actual);
ss2=size(approx);
if any(ss1 ~= ss2)
  disp('error in multiangle: dim"s don"t match.')
  return;
end
N = ss1(1);
angles = zeros(N,1);
for i=1:N
    nn1 = norm(actual(i,:));
    nn2 = norm(approx(i,:));
    if ~(nn1==0 || nn2==0)
        tt = actual(i,:)*approx(i,:)'/(norm(actual(i,:))*norm(approx(i,:)));
        if tt<(-1-1e-4) || tt>(+1+1e-4)
            angles(i) = NaN;
        else
            tt = mycliplims(tt,-1.0,+1.0);
            angles(i) = 180/pi*acos( tt);
        end
    end
end

