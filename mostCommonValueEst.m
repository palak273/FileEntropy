function entEst = mostCommonValueEst(data,len,alph_size)
ZALPHA = 2.5758293035489008;
for i=1:alph_size
    counts(i) = 0;
end
for i=1:len
    counts(data(i)+1) = counts(data(i)+1)+1;
end
mode = 0;
for i=1:alph_size
    if(counts(i)>mode) 
        mode = counts(i);
    end
end
pmax = mode/len;
ubound = min(1.0,pmax + ZALPHA*sqrt(pmax*(1.0-pmax)/(len-1.0)));
entEst = -log2(ubound);
end