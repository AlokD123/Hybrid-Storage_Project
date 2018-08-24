%Get:
%   1) minimum q-value for each state-action, used to determine optimal control
%   2) number of optimal controls for each state, based on NUMBER OF
%   SIMILAR Q-VALUES for all actions in that state ... NumOptCtrls_1
%   3) number of optimal controls for each state, based on POLICY MATRIX
%   (second-check)... NumOptCtrls_2

NumOptCtrls_1=[];
NumOptCtrls_2=[];
MaxQValues=[];

l=1;
for i=1:N1
    for j=1:N2
        for k=1:size(feasStates,3)
            if(feasStates(i,j,k)==1)
               MinQValues=[MaxQValues;i,j,k,min(qValsMtx(i,j,k,:))];
               NumOptCtrls_1=[NumOptCtrls_1;i,j,k,nnz(abs(qValsMtx(i,j,k,:)-min(qValsMtx(i,j,k,:)))<epsilon5)];
               l=l+1;
               NumOptCtrls_2=[NumOptCtrls_2;i,j,k,nnz(fullPolicyMtx(i,j,k,:))];
            end
        end
    end
end