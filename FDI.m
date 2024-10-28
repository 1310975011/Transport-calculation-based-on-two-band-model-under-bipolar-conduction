function y =FDI(i,n)
    fun = @(x)x.^i./(1+exp(x-n));
    y = integral(fun,0,inf,'ArrayValued',true);
end