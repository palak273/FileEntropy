function entEst = predictionEstimate(C,N,max_run_len,k) 
ZALPHA = 2.5758293035489008;
p_local= -1.0;
curMax = 1.0 / k;
p_global = C/N;
if(p_global > 0) 
    p_globalPrime = min(1.0, p_global + ZALPHA*sqrt((p_global*(1.0-p_global))/(N-1.0)));
else
    p_globalPrime = 1.0 - power(0.01, 1.0/N);
end
curMax = max(curMax, p_globalPrime);
if((curMax < 1.0) && (prediction_estimate_function(curMax, max_run_len+1, N) > log(0.99))) 
    p_local = calc_p_local(max_run_len, N, curMax);
    curMax = max(curMax, p_local);
end
entEst = -log2(curMax);
end
