function [prob_nextE] = ProbDistr_v2(nnzProb_nextE,currState_indL,nextState_indLs_Vect,resL_Mult)
%ProbDistr Transforms the input probabilities to a custom distribution
%CUSTOM DISTRIBUTION: BINOMIAL
%   Input: length of vector of uniform probabilities of next states (next demands);
%           indL in single current state;
%           vector of all possible indLs in next E-state;
%           demand resolution multiplicative factor (natural number)
%   Output: probabilities, with higher weight on staying near same
%           state (demand)

global epsilon3;

%Transform #of next loads input by multiplicative factor to achieve higher
%resolution of demand (during simulation, with continuous loads)
nnzProb_nextE=resL_Mult*(nnzProb_nextE-1)+1;


%Set number of "samples" controls number of possible next state loads (finite, non-negative)
%Binomial PDF parameter n
n=nnzProb_nextE-1; %One less because lowest binomial random variable is 0 instead of 1

%Small term to avoid all-none probability cases when n!=0 in binomial distribution
epsP=(1/n)/2; %Change p by a term less than the resolution of p steps (1/n), so not making other load more probable

%Find number of elements to shift in weights (below) to allow for HIGHEST WEIGHT ON RETURNING TO SAME STATE (LOAD)
%i.e. highest weight on MAIN DIAGONAL
d=find(abs(currState_indL-nextState_indLs_Vect)<epsilon3);
%If current load not possible in next state (under ANY control), find closest feasible one
diffL=1; %difference between current load and closest feasible one next
while isempty(d)
    d=find(abs(currState_indL-nextState_indLs_Vect)<diffL);
    diffL=diffL+1; %Step by discretization in load %<--------------------------------------- TO DO: make <1 if deltaL discretization is <1 (almost never)
end

%Transform d in the same way as nnzProb_nextE
d=resL_Mult*(d-1)+1;

%Create vector containing BINOMIAL PDF PROBABILITIES <---------------------------------------------- DISTRIBUTION W/ BINOMIAL PDF!!!!
%Set binomial distribution parameter "p", which controls probability distribution over all "n" feasible loads
if n==0 %If only one next state, value of "p" doesn't matter
    p=0.5; %Arbitrary
    
%If all or no shift, only have one possible load under this scheme. So...
elseif (d-1)/n==1   %If all shift, make smaller loads possible
    p=(d-1)/n-epsP;  %Decrease by small term
elseif (d-1)/n==0
    p=(d-1)/n+epsP;  %In opposite case, increase by small term
else
    p=(d-1)/n; %Set to be fraction of "n" so that shift factor towards larger loads shifts distribution
end
binomPDF_probs=binopdf(0:n,n,p);

%Return custom probabilities
prob_nextE=binomPDF_probs;

end

