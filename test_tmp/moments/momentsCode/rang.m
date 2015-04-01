function therange  = rang(mat)
s = size(mat);
mx = squeeze(max(mat));
themax = max(mx);
mn = squeeze(min(mat));
themin = min(mn);
therange = [ themin themax ] ;
