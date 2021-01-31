function a = G(z,d,num_blocks)
Ai=0.0; Ai_comp=0.0;
firstSum=0.0; firstSum_comp=0.0;
v = num_blocks - d;
Bterm = (1.0-z);
Bi = Bterm;
for i=2:d
    [Ai, Ai_comp] = kahan_add(Ai, Ai_comp, log2(i)*Bi);
    Bi = Bi*Bterm;
end
Ad1 = Ai; underflowTruncate = false;
for i=d+1:num_blocks-1
    ai = log2(i)*Bi;
    [Ai, Ai_comp] = kahan_add(Ai, Ai_comp, ai);
    aiScaled = (num_blocks-i) * ai;
    if(aiScaled > 0.0) 
        [firstSum, firstSum_comp] = kahan_add(firstSum,firstSum_comp,aiScaled);
    else
        underflowTruncate = true;
        break;
    end
    Bi = Bi*Bterm;
end
[firstSum, firstSum_comp] = kahan_add(firstSum,firstSum_comp,(num_blocks-d)*Ad1);
if(~underflowTruncate) 
    ai = log2(num_blocks)*Bi;
    [Ai, Ai_comp] = kahan_add(Ai,Ai_comp,ai);
end
a = 1/v * z*(z*firstSum + (Ai - Ad1));
end