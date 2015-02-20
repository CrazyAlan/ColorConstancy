function tm = trimean(x)
tm = ( median(x)*2+quantile(x,0.25)+quantile(x,0.75) )/4;