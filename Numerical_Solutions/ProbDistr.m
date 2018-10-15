function [prob_nextE] = ProbDistr(unifProbs_nxtState,nnzProb_nextE,currState_indL,nextState_indLs_Vect)
%ProbDistr Transforms the input probabilities to a custom distribution
%CUSTOM DISTRIBUTION: LINEAR PDF
%   Input: vector of uniform probabilities of next states (next demands);
%           length of said vector; indL in current state;
%           vector of all possible indLs in next E-state
%   Output: weighted probabilities, with higher weight on staying near same
%           state (demand)

global epsilon3

%Create convex combination vector containing WEIGHTS TO PRODUCT LINEAR PDF <---------------------------------------------- DISTRIBUTION W/ LINEAR PDF!!!!
linPDF_weights=(nnzProb_nextE:-1:1);
%Find number of elements to shift in linPDF_weights to allow for HIGHEST WEIGHT ON RETURNING TO SAME STATE (LOAD)
%i.e. highest weight on MAIN DIAGONAL
d=find(abs(currState_indL-nextState_indLs_Vect)<epsilon3);
%If current load not possible in next state (under ANY control), find closest feasible one
diffL=1; %difference between current load and closest feasible one next
while isempty(d)
    d=find(abs(currState_indL-nextState_indLs_Vect)<diffL);
    diffL=diffL+1; %Step by discretization in load %<--------------------------------------- TO DO: make <1 for interpolation load
end
%Shift vector of weights to right (highest to diagonal)
linPDF_weights=circshift(linPDF_weights,d-1);
%ATTEMPT TO MAKE PROBABILITIES INCREASE CLOSEST TO MAIN DIAGONAL OF
%P-MATRIX: shift all but main diagonal element back to left
if d~=1
    offDiag_weights=linPDF_weights(1:end~=d); %Get vector of off-diagonal elements
    linPDF_weights(1:end~=d)=circshift(offDiag_weights,-1); %Shift left and re-assign
end
%Transform original probabilities
prob_nextE=unifProbs_nxtState.*linPDF_weights'/sum(unifProbs_nxtState.*linPDF_weights');

end

