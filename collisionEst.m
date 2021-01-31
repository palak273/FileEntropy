function entEst = collisionEst(data,len)
ZALPHA = 2.5758293035489008;
i = 0; v = 0; s = 0.0;
while (i<len-1)
    if (data(i+1)==data(i+2)) 
        t_v = 2; %00 or 11
    elseif (i<len-2) 
        t_v = 3; %101, 011, 100, or 101
    else
        break;
    end
    v = v+1;
    s = s+t_v*t_v;
    i = i+t_v;
end
%X is mean of t_v's, s is sample stdev, where
%s^2 = (sum(t_v^2) - sum(t_v)^2/v) / (v-1)
X = i/v;
s = sqrt((s-(i*X))/(v-1));
X = X-ZALPHA * s/sqrt(v);
%2 is the smallest meaninful value here.
if(X<2.0)
    X = 2.0;
end
if(X<2.5)
    p = 0.5 + sqrt(1.25 - 0.5 * X);
    entEst = -log2(p);
else
    p = 0.5;
    entEst = 1.0;
end
end