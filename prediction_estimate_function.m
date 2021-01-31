function a = prediction_estimate_function(p,r,N)
LDBL_EPSILON = 1.0000e-09;
xlast=0.0;
q = 1.0-p;
x = 1.0;
i=0;
while (i <= 65) && ((x - xlast) > (LDBL_EPSILON*x)) 
    xlast = x;
    x = 1.0 + q*power(p, r)*power(x, r+1.0);
    i=i+1;
end
a = ((log(1.0-p*x) - log((r+1.0-r*x)*q) - (N+1.0)*log(x)));
end

