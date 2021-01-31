function [sum,comp] = kahan_add(sum,comp,in)
	y = in - comp;
	t = sum + y;
	comp = (t - sum) - y;
	sum = t;
end