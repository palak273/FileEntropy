clc
clear all
dp.word_size = 0;
initial_entropy = true;
all_bits = true;
file_path = 'test.txt';
MIN_SIZE = 1000000;
%%%%read_file%%%%
file = fopen(file_path);
if(file==-1)
    disp("Error: could not open '%s'\n", file_path);
end
rc = fseek(file, 0,'eof');
if(rc==-1)
    printf("Error: fseek failed\n");
    fclose(file);
end
dp.len = ftell(file);
if(dp.len<0)
    disp("Error: ftell failed\n");
    fclose(file);
end
frewind(file);
if(dp.len == 0)
    disp("Error: '%s' is empty\n", file_path);
    fclose(file);
end
dp.symbols = fread(file,dp.len);
fclose(file);
if(dp.word_size == 0) 
    datamask = 0; curbit = 128;
    for i=1:dp.len
        datamask = bitor(datamask,dp.symbols(i));
    end
    i=8;
    while (i>0) && (bitand(datamask,curbit) == 0)
        curbit = bitshift(curbit,-1);
        i = i-1;
    end
    dp.word_size = i;
else
    datamask = 0; curbit = 128;
    for i=1:dp.len
        datamask = bitor(datamask, dp.symbols(i));
    end
    i=8;
    while (i>0) && ((datamask & curbit) == 0)
        curbit = bitshift(curbit,-1);
        i = i-1;
    end
    if i<dp.word_size
        disp("Warning: Symbols appear to be narrower than described.\n");
    elseif i>dp.word_size
        disp("Incorrect bit width specification: Data does not fit within described bit width.\n");
        dp.symbols = NULL;
        dp.rawsymbols = NULL;
    end
end
dp.rawsymbols = typecast(dp.symbols,'single');
%memcpy(dp->rawsymbols, dp->symbols, sizeof(byte)* dp->len);%
dp.maxsymbol = 0;
max_symbols = bitshift(1,dp.word_size);
%int symbol_map_down_table[max_symbols];
dp.alph_size = 0;
%memset(symbol_map_down_table, 0, max_symbols*sizeof(int));
symbol_map_down_table = zeros([1 max_symbols]);
mask = max_symbols-1;
for i=1:dp.len 
    dp.symbols(i) = bitand(dp.symbols(i),mask);
    if (dp.symbols(i)>dp.maxsymbol) 
        dp.maxsymbol = dp.symbols(i);
    end
    if (symbol_map_down_table(dp.symbols(i)) == 0) 
        symbol_map_down_table(dp.symbols(i)) = 1;
    end
end
for i=1:max_symbols
    if (symbol_map_down_table(i) ~= 0) 
        symbol_map_down_table(i) = dp.alph_size;
        dp.alph_size = dp.alph_size+1;
    end
end
dp.blen = dp.len * dp.word_size;
if(dp.word_size == 1) 
    dp.bsymbols = dp.symbols;
else
    for i=1:dp.len
        for j=1:dp.word_size
            dp.bsymbols(i*dp.word_size+j) = bitand(bitshift(dp.symbols(i),-(dp.word_size-1-j)),1);
        end
    end
end
if(dp.alph_size < dp.maxsymbol + 1)
    for i=1:dp.len 
        dp.symbols(i) = symbol_map_down_table(dp.symbols(i));
    end
end
%%%%%%%%%%%%%%%%%

if(dp.alph_size <= 1)
    disp("Symbol alphabet consists of 1 symbol. No entropy awarded...\n");
end
if(~all_bits && (dp.blen>MIN_SIZE))
    dp.blen = MIN_SIZE;
end
if(((dp.alph_size > 2) || ~initial_entropy)) 
    disp("Number of Binary Symbols: %ld\n");
    disp(dp.blen);
end
if(dp.len<MIN_SIZE) 
    disp("\n*** Warning: data contains less than %d samples ***\n\n");
end
if (dp.alph_size < bitshift(1,dp.word_size)) 
    disp("\nSymbols have been translated.\n");
end
disp("Calculating Entropy");
H_original = dp.word_size; H_bitstring = 1.0;
if(((dp.alph_size > 2) || ~initial_entropy)) 
    ret_min_entropy = mostCommonValueEst(dp.bsymbols, dp.blen, 2);
    H_bitstring = min(ret_min_entropy, H_bitstring);
    disp("Most Common Value Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_bitstring = "+ H_bitstring);
end
if (initial_entropy) 
    ret_min_entropy = mostCommonValueEst(dp.symbols, dp.len, dp.alph_size);
    H_original = min(ret_min_entropy, H_original);
    disp("Most Common Value Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_original = "+ H_original);
end
if(((dp.alph_size > 2) || ~initial_entropy)) 
    ret_min_entropy = collisionEst(dp.bsymbols, dp.blen);
    H_bitstring = min(ret_min_entropy, H_bitstring);
    disp("Collision Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_bitstring = "+ H_bitstring);
end
if(initial_entropy && (dp.alph_size == 2))
    ret_min_entropy = collisionEst(dp.symbols, dp.len);
    H_original = min(ret_min_entropy, H_original);
    disp("Collision Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_original = "+ H_original);
end
if(((dp.alph_size > 2) || ~initial_entropy)) 
    ret_min_entropy = markovEst(dp.bsymbols, dp.blen);
    H_bitstring = min(ret_min_entropy, H_bitstring);
    disp("Markov Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_bitstring = "+ H_bitstring);
end
if(initial_entropy && (dp.alph_size == 2)) 
    ret_min_entropy = markovEst(dp.symbols, dp.len);
    H_original = min(ret_min_entropy, H_original);
    disp("Markov Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_original = "+ H_original);
end
if(((dp.alph_size > 2) || ~initial_entropy)) 
    ret_min_entropy = compressionEst(dp.bsymbols,dp.blen);
    if(ret_min_entropy >= 0)
        H_bitstring = min(ret_min_entropy, H_bitstring);
    end
    disp("Compression Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_bitstring = "+ H_bitstring);
end
if(initial_entropy && (dp.alph_size == 2)) 
    ret_min_entropy = compressionEst(dp.symbols,dp.len);
    H_original = min(ret_min_entropy, H_original);
    disp("Compression Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_original = "+ H_original);
end
if(((dp.alph_size > 2) || ~initial_entropy)) 
    ret_min_entropy = multiMCWTest(dp.bsymbols, dp.blen, 2);
    if (ret_min_entropy >= 0)
        H_bitstring = min(ret_min_entropy, H_bitstring);
    end
    disp("MCW Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_bitstring = "+ H_bitstring);
end
if(initial_entropy)
    ret_min_entropy = multiMCWTest(dp.symbols, dp.len, dp.alph_size);
    if(ret_min_entropy >= 0)
        H_original = min(ret_min_entropy, H_original);
    end
    disp("MCW Estimate");
    disp("ret_min_entropy = "+ret_min_entropy+" H_original = "+ H_original);
end