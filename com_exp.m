function a = com_exp(p,alph_size,d,num_blocks)
q = (1.0-p)/(alph_size-1.0);
a = G(p,d,num_blocks)+(alph_size-1.0)*G(q,d,num_blocks);
end