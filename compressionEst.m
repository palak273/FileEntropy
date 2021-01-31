function entEst = compressionEst(data,len)
ZALPHA = 2.5758293035489008;
ITERMAX = 1076;
DBL_INFINITY = inf;
RELEPSILON = 1.0000e-09;
ABSEPSILON = 1.0000e-37;
b = 6;
alph_size = bitshift(1,b); 
X=0.0; X_comp=0.0;
sigma=0.0; sigma_comp=0.0;
d = 4;
num_blocks = len/b;
if(num_blocks <= d)
    disp("\t*** Warning: not enough samples to run compression test (need more than %d) ***\n");
    entEst = -1.0;
    return;
end
for i=0:alph_size-1
    dict(i+1) = 0;
end
for i=0:d-1
    block = 0;
    for j=0:b-1
        block = bitor(block,bitshift(bitand(data(i*b+j+1),1),(b-j-1)));
		dict(block+1) = i+1;
    end
end
v = num_blocks - d;
for i=d:num_blocks-1
    block = 0;
    for j=0:b-1
        block = bitor(block,bitshift(bitand(data(i*b+j+1),1),(b-j-1)));
    end
    [X,X_comp] = kahan_add(X, X_comp, log2(i+1-dict(block+1)));
    [sigma,sigma_comp] = kahan_add(sigma, sigma_comp, log2(i+1-dict(block+1))*log2(i+1-dict(block+1)));
    dict(block+1) = i+1;
end
%compute mean and stdev
X = X/v;
sigma = 0.5907 * sqrt(sigma/(v-1.0) - X*X);
%binary search for p
X = X-ZALPHA * sigma/sqrt(v);
if(com_exp(1.0/alph_size,alph_size,d,num_blocks)>X)
    ldomain = 1.0/alph_size; hdomain = 1.0;
    lbound = ldomain; hbound = hdomain;
    lvalue = DBL_INFINITY; hvalue = -DBL_INFINITY;
    p = (lbound+hbound)/2.0;
    pVal = com_exp(p,alph_size,d,num_blocks);
    for j=0:ITERMAX-1
        if(relEpsilonEqual(pVal, X, ABSEPSILON, RELEPSILON, 4)) 
            break;
        end
        if(X < pVal) 
            lbound = p; lvalue = pVal;
        else
            hbound = p; hvalue = pVal;
        end
        if(lbound >= hbound)
            p = min(max(lbound, hbound),hdomain);
            break;
        end
        if(~(INCLOSEDINTERVAL(lbound, ldomain, hdomain) && INCLOSEDINTERVAL(hbound,  ldomain, hdomain)))
            p = ldomain;
            break;
        end
        if(~INCLOSEDINTERVAL(X, lvalue, hvalue)) 
            p = ldomain;
            break;
        end
        lastP = p;
        p = (lbound + hbound) / 2.0;
        if(~INOPENINTERVAL(p,lbound,hbound)) 
            p = hbound;
            break;
        end
        if(lastP == p) 
            p = hbound;
            break;
        end
        pVal = com_exp(p, alph_size, d, num_blocks);
        if(~INCLOSEDINTERVAL(pVal, lvalue, hvalue)) 
            p = hbound;
            break;
        end
    end
else
    p = -1.0;
end
if(p>1.0/alph_size) 
    entEst = -log2(p)/b;
else
    p = 1.0/alph_size;
    entEst = 1.0;
end
end