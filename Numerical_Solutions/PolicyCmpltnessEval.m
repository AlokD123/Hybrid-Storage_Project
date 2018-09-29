%Count number of optimal controls for each state

NumOptCtrls=[];
epsLM=1e-2;

for i=1:N1
    for j=1:N2
        for k=1:size(feasStates,3)
            if(feasStates(i,j,k)==1)
                %Count # of optimal controls in each state as those where
                %lagrange multipliers are >
               NumOptCtrls=[NumOptCtrls;i,j,k,nnz(fullPolicyMtx(i,j,k,:)>epsLM)];
               %NumOptCtrls=[NumOptCtrls;i,j,k,max((fullPolicyMtx(i,j,k,:)))];
            end
        end
    end
end