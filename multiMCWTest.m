function entEst = multiMCWTest(data,len,alph_size)
NUM_WINS = 4;
% W = [63, 255, 1023, 4095];
W = [3 5 7 9]; 
scoreboard = zeros(1,NUM_WINS);
max_cnts = zeros(1,NUM_WINS);
if(len<W(NUM_WINS)+1)	
    disp("\t*** Warning: not enough samples to run multiMCW test (need more than %d) ***\n");
    entEst = -1.0;
    return;
end
N = len-W(1);
winner = 0;
C = 0;
run_len = 0;
max_run_len = 0;
for i=0:NUM_WINS-1
    for j=0:alph_size-1
        win_cnts(i+1,j+1) = 0;
        win_poses(i+1,j+1) = 0;
    end
end
for i=0:W(NUM_WINS)-1
    for j=0:NUM_WINS-1
        if i<W(j+1)
            win_cnts(j+1,data(i+1)+1) = win_cnts(j+1,data(i+1)+1)+1;
            if max_cnts(j+1) <= win_cnts(j+1,data(i+1)+1)
                max_cnts(j+1) = win_cnts(j+1,data(i+1)+1);
                frequent(j+1) = data(i+1);
            end
            win_cnts(j+1,data(i+1)+1) = i;
        end
    end
end
for i=W(1):len-1
    if(frequent(winner+1) == data(i+1))
        C = C+1;
        run_len = run_len+1;
        if(run_len > max_run_len) 
            max_run_len = run_len;
        end
    else
        run_len = 0;
    end
    for j=0:NUM_WINS-1
        if((i >= W(j+1)) && (frequent(j+1) == data(i+1)))
            scoreboard(j+1) = scoreboard(j+1)+1;
            if(scoreboard(j+1) >= scoreboard(winner+1)) 
                winner = j;
            end
        end
    end
    for j=0:NUM_WINS-1
        if(i >= W(j+1))
            win_cnts(j+1,data(i-W(j+1)+1)+1) = win_cnts(j+1,data(i-W(j+1)+1)+1)-1;
            win_cnts(j+1,data(i+1)+1) = win_cnts(j+1,data(i+1)+1)+1;
            win_poses(j+1,data(i+1)+1) = i+1;
            if((data(i-W(j+1)+1) ~= frequent(j+1)) && (max_cnts(j+1) <= win_cnts(j+1,data(i+1)+1)))
                max_cnts(j+1) = win_cnts(j+1,data(i+1)+1);
                frequent(j+1) = data(i+1);
            elseif (data(i-W(j+1)+1) == frequent(j+1))
                max_cnts(j+1) = max_cnts(j+1)-1;
                max_pos = i-W(j+1);
                for k=0:alph_size-1
                    if((max_cnts(j+1) < win_cnts(j+1,k+1)) || ((max_cnts(j+1) == win_cnts(j+1,k+1)) && (max_pos <= win_poses(j+1,k+1))))
                        max_cnts(j+1) = win_cnts(j+1,k+1);
                        frequent(j+1) = k;
                        max_pos = win_poses(j+1,k+1);
                    end
                end
            end
        end
    end
end
entEst = predictionEstimate(C, N, max_run_len, alph_size);
end