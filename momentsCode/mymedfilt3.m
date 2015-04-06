function B = mymedfilt3(A,siz)

B(:,:,1) = medfilt2(squeeze(A(:,:,1)),siz);
B(:,:,2) = medfilt2(squeeze(A(:,:,2)),siz);
B(:,:,3) = medfilt2(squeeze(A(:,:,3)),siz);

