function p = calc_p_local(max_run_len,N,ldomain)
ITERMAX = 1076;
DBL_INFINITY = inf;
RELEPSILON = 1.0000e-09;
ABSEPSILON = 1.0000e-37;
log_alpha = log(0.99);
hdomain = 1.0;
lbound = ldomain; hbound = hdomain;
lvalue = DBL_INFINITY; hvalue = -DBL_INFINITY;
p = (lbound + hbound) / 2.0;
pVal = prediction_estimate_function(p, max_run_len+1, N);
for j=0:ITERMAX-1
    if(relEpsilonEqual(pVal, log_alpha, ABSEPSILON, RELEPSILON, 4)) 
        break;
    end
    if(log_alpha < pVal)
        lbound = p;
        lvalue = pVal;
    else
        hbound = p;
        hvalue = pVal;
    end
    if(lbound >= hbound) 
        p = min(max(lbound, hbound),hdomain);
        break;
    end
    if(~(INCLOSEDINTERVAL(lbound, ldomain, hdomain) && INCLOSEDINTERVAL(hbound,  ldomain, hdomain))) 
        p = hdomain;
        break;
    end
    if(~INCLOSEDINTERVAL(log_alpha, lvalue, hvalue))
        p = hdomain;
        break;
    end
    lastP = p;
    p = (lbound + hbound) / 2.0;
    if(~INOPENINTERVAL(p,  lbound, hbound)) 
        p = hbound;
        break;
    end
    if(lastP == p) 
        p = hbound;
        break;
    end
    pVal = prediction_estimate_function(p, max_run_len+1, N);
    if(~INCLOSEDINTERVAL(pVal, lvalue, hvalue)) 
        p = hbound;
        break;
    end
end
disp(p);
end
