% do_graham_moments_gehler_edges.m

% from: for_appletalkApr2012_canon1d.m

% from:  do_canon1d_diividebylight.m (from do_canon1d.m)
%
clear all;
%get_stuart_canon5D_filelist; % 482 images
%   NO: writesmallGehlerImages;
% makesmallGehlerImages;  % makes allcanon5Dsmall
%    readgehler;
load ('/home/xca64/remote/GitHub/colorP/dataSet/sfudataset/allImageMedFilter4.mat') % allcanon5Dsmall = zeros(482,183,275,3); % all portrait
allcanon5Dsmall = allImage;clear allImage;
[howmanycands, r,c, n3] = size(allcanon5Dsmall); % 482   183   275 3
% getGehlerLights; % gets alllightschrom

load('/home/xca64/remote/GitHub/colorP/dataSet/sfudataset/sfuIllum.mat');
alllights=sfuIllum; clear sfuIllum % 482 5DCimages

% what is grey? (for this camera):
alllightschrom3 = makechrom3vec(alllights);  %Normalize the illuminants ? 

clear D M C W % M=MomentsArray
% and solve D(diag 86x86) * M(86x9) * C(9x3)=W(whitepint chroms, 86x3)
%   for D and C:
dim=18
% M = zeros(howmanycands,dim);
allim = reshape(allcanon5Dsmall,howmanycands,r*c,3);
trainbot=1;
traintop=floor(2*howmanycands/3) % 321

mink_norm = 7;
sigma = 0.75;
diff_order = 1;

K = 3;
run_times = 5;
results = zeros(K*run_times,5);

%for imsize = 0.1:0.1:0.6
%    for sigma = 0.25:0.5:2.5
imsize = 0.3;
sigma = 2;

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
        
  %  H = fspecial('gaussian', [3 3], sigma) ; % default value for hsize is [3 3]; the default value for sigma is 0.5. 
 %   tmp = imfilter(reshape(im,r,c,3),H,'replicate');
 %   tmp = mymedfilt3(reshape(im,r,c,3),[5,5]);
  %   imsmooth = imresize(tmp,imsize);
%    [r_reshaped,c_reshaped, n3_reshaped] = size(imsmooth);
       imsmooth = reshape(im,r,c,3);
        [dx,dy] = makegradient(imsmooth);
    %   [dx, dy] = norm_derivative_separate_dxdy(imsmooth,sigma,1);
      %  imedges= sqrt( dx.^2+dy.^2 );
      %  imedges=reshape(imedges,r*c,3);
     %    dx=reshape(dx,r_reshaped*c_reshaped,3);
    %     dy=reshape(dy,r_reshaped*c_reshaped,3);
     dx=reshape(dx,r*c,3);
     dy=reshape(dy,r*c,3);
        M(i,:) = moments9_edges(dx,dy); % 65 x 18

    
       % M(i,:) = moments9(imedges);
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
        
%        imedges= sqrt( dx.^2+dy.^2 );
%        imedges=reshape(imedges,r*c,3);

   %    H = fspecial('gaussian', [3 3], sigma) ; % default value for hsize is [3 3]; the default value for sigma is 0.5. 
  %    tmp = imfilter(reshape(im,r,c,3),H,'replicate');
 %  tmp = mymedfilt3(reshape(im,r,c,3),[5,5]);
 %    imsmooth = imresize(tmp,imsize);
    %    [r_reshaped,c_reshaped, n3_reshaped] = size(imsmooth);
        imsmooth = reshape(im,r,c,3);
        [dx,dy] = makegradient(imsmooth);
    %   [dx, dy] = norm_derivative_separate_dxdy(imsmooth,sigma,1);
      %  imedges= sqrt( dx.^2+dy.^2 );
      %  imedges=reshape(imedges,r*c,3);
      %   dx=reshape(dx,r_reshaped*c_reshaped,3);
       %  dy=reshape(dy,r_reshaped*c_reshaped,3);
       dx=reshape(dx,r*c,3);
       dy=reshape(dy,r*c,3);
        
        anM = moments9_edges(dx,dy); % 65 x 18
       % anM = moments9(imedges);
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
%mea(int8(imsize*10),int8((sigma*4)/2),:) = mean(results)
%std1(int8(imsize*10),int8((sigma*4)/2),:) = std(results)
%    end
%end
% 
% %% PRIOR ART:
% allpriorerrs=zeros(howmanycands,4);
% for iii=1:howmanycands % iii=22
%     disp(iii);
%     im = squeeze(allcanon5Dsmall(iii,:,:,:));
%     im(im==0) = 1/256;
%     [r,c,n3]=size(im); % 256 170 3
%     correctii = iii;
%              maxrgb = zeros(3,1);
%             greyrgb = zeros(3,1);
%             minkrgb = zeros(3,1); pnorm=6;
%             greyedgergb = zeros(3,1);
%              [dx,dy] = makegradient(im);
%             for k=1:3
%                 maxrgb(k)=max(max(im(:,:,k)));
%                 greyrgb(k) = mean(mean(im(:,:,k)));
%                 minkrgb(k) = (1/(r*c)^(1/pnorm)) * (sum(sum( (im(:,:,k).^pnorm))))^(1/pnorm);
%                 greyedgergb(k)=mean(mean( sqrt( dx(:,:,k).^2+dy(:,:,k).^2 ) ));
%             end
%             maxrgbchrom = maxrgb/sum(maxrgb); % 0.4010    0.4101    0.1889
%             greyrgbchrom = greyrgb/sum(greyrgb); % 0.3734    0.4184    0.2082
%             minkrgbchrom = minkrgb/sum(minkrgb); % 0.3734    0.4184    0.2082
%             greyedgergbchrom = greyedgergb/sum(greyedgergb);
%             allpriorerrs(iii,:)=( multiangle(repmat(alllightschrom3(correctii,:),[4 1]),...
%                 [maxrgbchrom';greyrgbchrom';minkrgbchrom';greyedgergbchrom']) )';%
% end; % for iii
% % RESULTS:
% [mean(allpriorerrs(:,4)), median(allpriorerrs(:,4)),...
%     trimean(allpriorerrs(:,4)), min(allpriorerrs(:,4)), ...
%     quantile(allpriorerrs(:,4),0.95)]
% % 4.7835    3.8740    4.1539    0.1092   11.1252


