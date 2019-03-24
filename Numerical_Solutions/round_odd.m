function S = round_odd(S)
% round to nearest odd integer.
idx = mod(S,2)<1;
S = floor(S);
S(idx) = S(idx)+1;