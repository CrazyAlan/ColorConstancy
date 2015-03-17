
% from: for_appletalkApr2012_canon1d.m

% from:  do_canon1d_diividebylight.m (from do_canon1d.m)
%
clear all;
%get_stuart_canon5D_filelist; % 482 images
%   NO: writesmallGehlerImages;
% makesmallGehlerImages;  % makes allcanon5Dsmall
%    readgehler;
load ('../dataSet/sfudataset/allImageOriginal.mat') % allcanon5Dsmall = zeros(482,183,275,3); % all portrait
allcanon5Dsmall = allImage;clear allImage;
[howmanycands, r,c, n3] = size(allcanon5Dsmall); % 482   183   275 3
% getGehlerLights; % gets alllightschrom

load('../dataSet/sfudataset/sfuIllum.mat');
alllights=sfuIllum; clear sfuIllum % 482 5DCimages


% what is grey? (for this camera):
alllightschrom3 = makechrom3vec(alllights);  %Normalize the illuminants ? 

clear D M C W % M=MomentsArray
% and solve D(diag 86x86) * M(86x9) * C(9x3)=W(whitepint chroms, 86x3)
%   for D and C:
dim=19
% M = zeros(howmanycands,dim);
allim = reshape(allcanon5Dsmall,howmanycands,r*c,3);
trainbot=1;
traintop=floor(2*howmanycands/3) % 321

K = 3;
run_times = 10;
results = zeros(K*run_times,5);

for t=1:run_times
    % Generate cross-validation indices
    cv_folds = crossvalind('Kfold', howmanycands, K); 

for k=1:K
    
    test_idx = (cv_folds == k);
    train_idx = ~test_idx;
    
    train_data = allim(train_idx,:,:);
    test_data = allim(test_idx,:,:);
    M = [];
    
    for i=1:size(train_data,1)
        im=squeeze(train_data(i,:,:));
        lum = sum(im,2);
        bright = lum>quantile(lum,0.95);
        
        [dx,dy] = makegradient(reshape(im,r,c,3));
        imedges= sqrt( dx.^2+dy.^2 );
        imedges=reshape(imedges,r*c,3);
    
        M(i,:) = moments19(imedges);
    end % for i
    W = alllightschrom3(train_idx,:);
    D=eye(size(M,1)); % init
    C=zeros(size(M,1),3);
    MaxItn=100;
    for i=1:MaxItn
        oldD=D;
        oldC=C;
        %
        C=pinv(D*M)*W;
        %C=( mytikhonov_fixedAmountofregularization( (D*M)' ,W' ,1e-9) )' ;
        for j=1:size(M,1)
            % LS:
            mc = M(j,:)*C; % 1x3
            D(j,j)=pinv( mc*mc') *  (mc*W(j,:)');
        end % for j
    end % for i

    % Ok, now test:
    chromout=zeros(size(test_data,1),3); % omit trainings later.
    chrom_truth = alllightschrom3(test_idx,:);
    allangerrs = zeros(size(test_data,1),1);
    % clear M
    for i=1:size(chrom_truth,1)
        im=squeeze(test_data(i,:,:));
        lum = sum(im,2);
        bright = lum>quantile(lum,0.95);
        
        [dx,dy] = makegradient(reshape(im,r,c,3));
        imedges= sqrt( dx.^2+dy.^2 );
        imedges=reshape(imedges,r*c,3);

        anM = moments19(imedges);
        chromout(i,:) = makechrom3vec( anM*C );
        allangerrs(i) =  multiangle(chrom_truth(i,:),chromout(i,:)); % degrees
    end % for i
   
   results((t-1)*K+k,:) =  [mean(allangerrs), median(allangerrs), trimean(allangerrs), ...
        min(allangerrs), quantile(allangerrs,0.95)];
    results((t-1)*K+k,:)
    % 3.2526    2.4019    2.5301    0.0396    8.5540
end

end
mea = mean(results)
std1 = std(results)