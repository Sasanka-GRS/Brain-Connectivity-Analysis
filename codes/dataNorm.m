function y = dataNorm(x)
    l = length(x);
    mean = sum(x)/l;
    var = x'*x/l - mean*mean;
    y = (x-mean)/var;
end