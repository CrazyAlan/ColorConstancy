clear all;

load /home/xca64/remote/GitHub/colorP/dataSet/sfudataset/allImageFourthResolutionNoFilter.mat
load /home/xca64/remote/GitHub/colorP/dataSet/sfudataset/sfuIllum.mat

[howmanycands, r,c, n3] = size(allImage); % 482   183   275 3
alllightschrom3 = makechrom3vec(sfuIllum);


allpriorerrs=zeros(howmanycands,4);
for iii=1:howmanycands % iii=22
    disp(iii);
    im = squeeze(allImage(iii,:,:,:));
    im(im==0) = 1/256;
    [r,c,n3]=size(im); % 256 170 3
    correctii = iii;
         maxrgb = zeros(3,1);
         greyrgb = zeros(3,1);
         minkrgb = zeros(3,1); pnorm=7;
         greyedgergb = zeros(3,1);
         [dx,dy] = makegradient(im);
            for k=1:3
                maxrgb(k)=max(max(im(:,:,k)));
                greyrgb(k) = mean(mean(im(:,:,k)));
                minkrgb(k) = (1/(r*c)^(1/pnorm)) * (sum(sum( (im(:,:,k).^pnorm))))^(1/pnorm);
%     sigma=5; % 5?
%     % imsmooth = imgaussfilt(reshape(im,r,c,3),sigma); % only in R2015 :(
%     H = fspecial('gaussian', [3 3], sigma) ; % default value for hsize is [3 3]; the default value for sigma is 0.5. 
%     imsmooth = imfilter(im(:,:,k),H,'replicate');
%     [dx,dy] = makegradient(imsmooth);
%     %imedges= sqrt( dx.^2+dy.^2 );
%     % imedges= ( abs(dx).^7+abs(dy).^7 ).^(1/7);
%     imedges= ( abs(dx).^2+abs(dy).^2 ).^(1/2);
%              %greyedgergb(k)=mean(mean( sqrt( dx(:,:,k).^2+dy(:,:,k).^2 ) ));
%               %  greyedgergb(k)=mean(imedges(:));
%              pnorm=7;
%             greyedgergb(k)= (1/(r*c)^(1/pnorm)) * (sum(sum( (imedges.^pnorm))))^(1/pnorm);
            end
       %   greyedgergb(k)= (1/(r*c)^(1/pnorm)) * (sum(sum( (imedges.^pnorm))))^(1/pnorm);
       
            mink_norm=7;    % any number between 1 and infinity
            sigma=1.5;        % sigma 
           diff_order=1;   % differentiation order (1 or 2)

           [wR,wG,wB,out4]=general_cc(im,diff_order,mink_norm,sigma);
%      [Rx,Gx,Bx] = norm_derivative(im,sigma,1);
%      input_data(:,:,1)=Rx;
%     input_data(:,:,2)=Gx;
%     input_data(:,:,3)=Bx;
%     input_data=abs(input_data);
%      kleur=power(input_data,mink_norm);
%      mask_im2 = ones(size(im,1),size(im,2));
%      mask_im2=set_border(mask_im2,sigma+1,0);
%      white_R = power(sum(sum(kleur(:,:,1).*mask_im2)),1/mink_norm);
%     white_G = power(sum(sum(kleur(:,:,2).*mask_im2)),1/mink_norm);
%     white_B = power(sum(sum(kleur(:,:,3).*mask_im2)),1/mink_norm);
% 
%     som=sqrt(white_R^2+white_G^2+white_B^2);
% 
%     white_R=white_R/som;
%     white_G=white_G/som;
%     white_B=white_B/som;
    
            greyedgergb = [wR;wG;wB];
          
            maxrgbchrom = maxrgb/sum(maxrgb); % 0.4010    0.4101    0.1889
            greyrgbchrom = greyrgb/sum(greyrgb); % 0.3734    0.4184    0.2082
            minkrgbchrom = minkrgb/sum(minkrgb); % 0.3734    0.4184    0.2082
            greyedgergbchrom = greyedgergb/sum(greyedgergb);
            allpriorerrs(iii,:)=( multiangle(repmat(alllightschrom3(correctii,:),[4 1]),...
                [maxrgbchrom';greyrgbchrom';minkrgbchrom';greyedgergbchrom']) )';%
end; % for iii
% RESULTS:
[mean(allpriorerrs(:,1)), median(allpriorerrs(:,1)),...
    trimean(allpriorerrs(:,1)), min(allpriorerrs(:,1)), ...
    quantile(allpriorerrs(:,1),0.95),...
    max(allpriorerrs(:,1))]
[mean(allpriorerrs(:,2)), median(allpriorerrs(:,2)),...
    trimean(allpriorerrs(:,2)), min(allpriorerrs(:,2)), ...
    quantile(allpriorerrs(:,2),0.95),...
    max(allpriorerrs(:,2))]
[mean(allpriorerrs(:,3)), median(allpriorerrs(:,3)),...
    trimean(allpriorerrs(:,3)), min(allpriorerrs(:,3)), ...
    quantile(allpriorerrs(:,3),0.95),...
    max(allpriorerrs(:,3))]
[mean(allpriorerrs(:,4)), median(allpriorerrs(:,4)),...
    trimean(allpriorerrs(:,4)), min(allpriorerrs(:,4)), ...
    quantile(allpriorerrs(:,4),0.95),...
    max(allpriorerrs(:,4))]