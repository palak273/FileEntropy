function entEst = markovEst(data,len)
C_0 = 0; C_00 = 0; C_10 = 0;
for i=0:len-2
    if(data(i+1)==0)
        C_0 = C_0+1;
        if(data(i+2) == 0) 
            C_00 = C_00+1;
        end
    elseif (data(i+2)==0)
        C_10 = C_10+1;
    end
end
%C_0 is now  the number of 0 bits from S[0] to S[len-2]
C_1 = len-1- C_0; %C_1 is the number of 1 bits from S[0] to S[len-2]
%Note that P_X1 = C_X1 / C_X = (C_X - C_X0)/C_X = 1.0 - C_X0/C_X = 1.0 - P_X0 
if(C_0 > 0) 
    P_00 = (C_00)/(C_0);
    P_01 = 1.0 - P_00;
else
    P_00 = 0.0;
    P_01 = 0.0;
end
if(C_1 > 0) 
    P_10 = (C_10)/(C_1);
    P_11 = 1.0 - P_10;
else
    P_10 = 0.0;
    P_11 = 0.0;
end
if (data(len)==0)
    C_0 = C_0+1;
	%C_0 is now  the number of 0 bits from S[0] to S[len-1]
end
P_0 = C_0/len;
P_1 = 1.0-P_0;
H_min = 128.0;
%In the next block, note that if P_0X > 0.0, then P_0 > 0.0
%and similarly if P_1X > 0.0, then P_1 > 0.0
%Sequence 00...0
if(P_00 > 0.0)
    tmp_min_entropy = -log2(P_0) - 127.0*log2(P_00); 
    if(tmp_min_entropy < H_min) 
        H_min = tmp_min_entropy;
    end
end
%Sequence 0101...01
if ((P_01 > 0.0) && (P_10 > 0.0))
    tmp_min_entropy = -log2(P_0) - 64.0*log2(P_01) - 63.0*log2(P_10);
    if(tmp_min_entropy < H_min) 
        H_min = tmp_min_entropy;
    end
end
%Sequence 011...1
if((P_01 > 0.0) && (P_11 > 0.0))
    tmp_min_entropy = -log2(P_0) - log2(P_01) - 126.0*log2(P_11);
    if(tmp_min_entropy < H_min)
        H_min = tmp_min_entropy;	
    end
end
%Sequence 100...0
if((P_10 > 0.0) && (P_00 > 0.0))
    tmp_min_entropy = -log2(P_1) - log2(P_10) - 126.0*log2(P_00);
    if(tmp_min_entropy < H_min) 
        H_min = tmp_min_entropy;
    end
end
%Sequence 1010...10
if((P_10 > 0.0) && (P_01 > 0.0))
    tmp_min_entropy = -log2(P_1) - 64.0*log2(P_10) - 63.0*log2(P_01);
    if(tmp_min_entropy < H_min)
        H_min = tmp_min_entropy;
    end
end
%Sequence 11...1
if(P_11 > 0.0)
    tmp_min_entropy = -log2(P_1) - 127.0*log2(P_11);
    if(tmp_min_entropy < H_min) 
        H_min = tmp_min_entropy;
    end
end
entEst = min(H_min/128.0,1.0);
end