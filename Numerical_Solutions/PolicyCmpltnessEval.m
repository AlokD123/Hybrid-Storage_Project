NumOptCtrls=[];

for i=1:N1
    for j=1:N2
        for k=1:size(feasStates,3)
            if(feasStates(i,j,k)==1)
               NumOptCtrls=[NumOptCtrls;i,j,k,nnz(fullPolicyMtx(i,j,k,:)>0)];
               %NumOptCtrls=[NumOptCtrls;i,j,k,max((fullPolicyMtx(i,j,k,:)))];
            end
        end
    end
end