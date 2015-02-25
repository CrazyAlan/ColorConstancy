% do_graham_moments_gehler.m

% from: for_appletalkApr2012_canon1d.m

% from:  do_canon1d_diividebylight.m (from do_canon1d.m)
%
get_stuart_canon5D_filelist; % 482 images
%   NO: writesmallGehlerImages;
% makesmallGehlerImages;  % makes allcanon5Dsmall
%    readgehler;
load 'allcanon5Dsmall.mat' allcanon5Dsmall % allcanon5Dsmall = zeros(482,183,275,3); % all portrait
[howmanycands, r,c, n3] = size(allcanon5Dsmall); % 482   183   275 3
% getGehlerLights; % gets alllightschrom
datapath = '/Users/crazyalan/GitHub/ColorConstancy/Gehler_Extras/'
load(strcat(datapath,'illuminants'));
alllights=illuminants; clear illuminants % 482 5DCimages

% what is grey? (for this camera):
alllightschrom3 = makechrom3vec(alllights);  %Normalize the illuminants ? 
% alllightschrom3mean = mean(alllightschrom3); %  0.2436    0.4276    0.3287
% alllightschrom3mean0 = zeros(howmanycands,3);
% for k=1:3
%     alllightschrom3mean0(:,k)=alllightschrom3(:,k)-alllightschrom3mean(k);
% end % now svd is 2D.

  %  alllightschromM= squeeze(makechromM_Nx3to3vec_exact(reshape(alllightschrom,howmanycands,1,3)));
    % Now log() is rank-2.

%{
% myplot3(log(alllightschromM(:,1:3))); % there is a line in there...
% axis equal
Pu=projector([1;1;1]);
Pu1 = eye(3)-Pu;
[U,D,V] = svd(Pu1);
V = V(:,1:2);
logalllightschromMPu = log(alllightschromM)*V;  % in-plane
[u,d,v]=svd(logalllightschromMPu,'econ');diag(d)
% quite line-like...
    myplot2(logalllightschromMPu); axis equal  %  :)
    [coeff,score,latent]=princomp(logalllightschromMPu); latent % 0.1891 0.0053
%}

%% Moments method:
%{
What Graham is actually using is: [and he uses RGB, or edges(RGB)]:
E(R) , where E() is the expectation (meaning the average);
E(G)
E(B)
( E(R^2) )^(1/2)
( E(G^2) )^(1/2)
( E(B^2) )^(1/2)
( E(R*G) )^(1/2)
( E(R*B) )^(1/2)
( E(G*B) )^(1/2) 
%}
clear D M C W % M=MomentsArray
% and solve D(diag 86x86) * M(86x9) * C(9x3)=W(whitepint chroms, 86x3)
%   for D and C:
dim=19
% M = zeros(howmanycands,dim);
allim = reshape(allcanon5Dsmall,howmanycands,r*c,3);
trainbot=1;
traintop=floor(2*howmanycands/3) % 321

A_idx = (1:floor(howmanycands/3));
B_idx = (ceil(howmanycands/3) : floor(2*howmanycands/3));
C_idx = (ceil(2*howmanycands/3) : howmanycands);

train_idx = [A_idx C_idx];
test_idx = [B_idx];

for i=1:size(train_idx,2)
    im=squeeze(allim(train_idx(i),:,:));
    lum = sum(im,2);
    bright = lum>quantile(lum,0.95);
    M(i,:) = moments19(im);
end % for i
W = alllightschrom3(train_idx,:);
D=eye(size(M,1)); % init
C=zeros(size(M,1),3);
MaxItn=100;
for j=1:MaxItn
    oldD=D;
    oldC=C;
    %
    C=pinv(D*M)*W;
    %C=( mytikhonov_fixedAmountofregularization( (D*M)' ,W' ,1e-9) )' ;
    for k=1:size(M,1)
        % LS:
        mc = M(k,:)*C; % 1x3
        %D(k,k)=pinv( mc*mc') *  (mc*W(k,:)'); % same as below:
        D(k,k)= W(k,:)*pinv(mc); % same as below:
    end % for k
end % for j

% Ok, now test:
testindices = 1:howmanycands;
testindices(trainbot:traintop)=[];
    ll=length(testindices);
chromout=zeros(howmanycands,3); % omit trainings later.
allangerrs = zeros(howmanycands,1);
% clear M
for i=test_idx
    im=squeeze(allim(i,:,:));
        lum = sum(im,2);
        bright = lum>quantile(lum,0.95);
            %bright=ones(r*c,1);
    %M(i,:) = moments9(im);
%     anM =[ moments9_geo(makechrom3vec(im(bright,:))) ...
%         moments9(im)];
    %M(i,:) = moments9_geo(makechrom3vec(im(:,:)));
anM = moments19(im);
    chromout(i,:) = makechrom3vec( anM*C );
    allangerrs(i) =  multiangle(alllightschrom3(i,:),chromout(i,:)); % degrees
end % for i
allangerrs(train_idx)=[];
hist(allangerrs)
% RESULTS:
[mean(allangerrs), median(allangerrs), trimean(allangerrs), ...
    min(allangerrs), quantile(allangerrs,0.95)]
% 3.2526    2.4019    2.5301    0.0396    8.5540


%% PRIOR ART:
allpriorerrs=zeros(howmanycands,4);
for iii=1:howmanycands % iii=22
    disp(iii);
    im = squeeze(allcanon5Dsmall(iii,:,:,:));
    im(im==0) = 1/256;
    [r,c,n3]=size(im); % 256 170 3
    correctii = iii;
             maxrgb = zeros(3,1);
            greyrgb = zeros(3,1);
            minkrgb = zeros(3,1); pnorm=6;
            greyedgergb = zeros(3,1);
             [dx,dy] = makegradient(im);
            for k=1:3
                maxrgb(k)=max(max(im(:,:,k)));
                greyrgb(k) = mean(mean(im(:,:,k)));
                minkrgb(k) = (1/(r*c)^(1/pnorm)) * (sum(sum( (im(:,:,k).^pnorm))))^(1/pnorm);
                greyedgergb(k)=mean(mean( sqrt( dx(:,:,k).^2+dy(:,:,k).^2 ) ));
            end
            maxrgbchrom = maxrgb/sum(maxrgb); % 0.4010    0.4101    0.1889
            greyrgbchrom = greyrgb/sum(greyrgb); % 0.3734    0.4184    0.2082
            minkrgbchrom = minkrgb/sum(minkrgb); % 0.3734    0.4184    0.2082
            greyedgergbchrom = greyedgergb/sum(greyedgergb);
            allpriorerrs(iii,:)=( multiangle(repmat(alllightschrom3(correctii,:),[4 1]),...
                [maxrgbchrom';greyrgbchrom';minkrgbchrom';greyedgergbchrom']) )';%
end; % for iii
% RESULTS:
[mean(allpriorerrs(:,4)), median(allpriorerrs(:,4)),...
    trimean(allpriorerrs(:,4)), min(allpriorerrs(:,4)), ...
    quantile(allpriorerrs(:,4),0.95)]
% 4.7835    3.8740    4.1539    0.1092   11.1252


