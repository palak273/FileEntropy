function ans = col_exp(p)
q = 1.0 - p;
ans = (p/(q*q))*(1.0 + 0.5*(1.0/p - 1.0/q))*F(q) - (p/q)*0.5*(1.0/p - 1.0/q);
end