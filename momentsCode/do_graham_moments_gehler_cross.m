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



clear D M C W % M=MomentsArray
% and solve D(diag 86x86) * M(86x9) * C(9x3)=W(whitepint chroms, 86x3)
%   for D and C:
dim=9
% M = zeros(howmanycands,dim);
allim = reshape(allcanon5Dsmall,howmanycands,r*c,3);
trainbot=1;
traintop=floor(2*howmanycands/3) % 321



for i=trainbot:traintop
    im=squeeze(allim(i,:,:));
        lum = sum(im,2);
        bright = lum>quantile(lum,0.95);
M(i,:) = moments9(im);
end % for i
W = alllightschrom3(trainbot:traintop,:);
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
for i=testindices
    im=squeeze(allim(i,:,:));
        lum = sum(im,2);
        bright = lum>quantile(lum,0.95);
anM = moments9(im);
    chromout(i,:) = makechrom3vec( anM*C );
    allangerrs(i) =  multiangle(alllightschrom3(i,:),chromout(i,:)); % degrees
end % for i
allangerrs(trainbot:traintop)=[];
hist(allangerrs)
% RESULTS:
[mean(allangerrs), median(allangerrs), trimean(allangerrs), ...
    min(allangerrs), quantile(allangerrs,0.95)]
% 3.2526    2.4019    2.5301    0.0396    8.5540



