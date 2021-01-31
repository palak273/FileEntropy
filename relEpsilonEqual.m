function che = relEpsilonEqual(A,B,maxAbsFactor,maxRelFactor,maxULP)
DBL_MIN = 1.0000e-37;
if(isnan(A) || isnan(B)) 
    che = false;
    return;
end
if(A==B)
    che = true;
    return;
end
if(isinf(A) || isinf(B)) 
    che = false;
    return ;
end
absA = abs(A);
absB = abs(B);
if(absA > absB)
    tmp = B;
    B=A;
    A=tmp;
    tmp = absB;
    absB = absA;
    absA = tmp;
end
diff=abs(B-A);
if((absA < DBL_MIN) || (diff < DBL_MIN) || isinf(diff) || (absB * maxRelFactor < DBL_MIN)) 
    che = (diff <= maxAbsFactor);
    return;
else
    if(diff <= absB * maxRelFactor) 
        che = true;
        return;
    end 
end
if(sign(A) ~= sign(B))
    che = false;
    return;
end
Aint = typecast(absA,'uint64');
Bint = typecast(absB,'uint64');
% memcpy(&Aint, &absA, sizeof(double));
% memcpy(&Bint, &absB, sizeof(double));
che = (Bint - Aint <= maxULP);
end