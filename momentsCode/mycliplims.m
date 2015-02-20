function imout=mycliplims(im,lo,hi)
imout=im;
imout(find(im<lo))=lo;
imout(find(im>hi))=hi;