function S = round2even(x) 
    if mod(x,2)<1 
        S = fix(x); 
    else 
        S =fix(x) + 1; 
    end
end