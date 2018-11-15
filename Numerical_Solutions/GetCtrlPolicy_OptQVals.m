%Create matrix of optimal q-values and use to obtain optimal control policy
%qValsMtx(E1,E2,L,:) has q-values of state (E1,E2,L) for all controls

global qValsMtx; global fullPolicyMtx;
qValsMtx=[]; fullPolicyMtx=[];

epsilon4=1e-3;

%N1=6; N2=5; M=10;
%P1=N1; P2=N2;

%Repeat for each decision (D1,D2 combination)
for p=1:p_max 
    %if (~isempty(Lmin{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...
        qSubVec=aug_Q((p-1)*N1*N2*M+1:p*N1*N2*M);
        %Format q-values sub-vector into E1xE2 matrices (M matrices, for each value of load)
            for i=0:M*N2:M*N2*(N1-1) 
                for j=0:M-1 
                  for ind=(1+i+j):M:(M*(N2-1)+i+(j+1)) 
                      if(mod(ind,M*N2)==0)
                        ind2=(M*N2-1-j)/M+1;
                      else
                        ind2=(mod(ind,M*N2)-1-j)/M+1;
                      end
                      qValsMtx(i/(M*N2)+1,ind2,j+1,p)=qSubVec(ind);
                  end
                end
            end
    %end
end

numOptCtrls=Inf; %TOTAL # of optimal ctrls. Initialize to infinity to guarantee 1 loop

while numOptCtrls > length(E_Ind_VectALL) %WHILE more than one optimal ctrl per state....
    epsilon4=epsilon4/10; %Decrease tolerance on minum Q-value detection UNTIL ONLY ONE OPT CTRL per STATE
    
    %NumOptCtrls_1=[];
    NumOptCtrls_2=[]; %Store number of optimal controls per state, based on FULL POLICY MTX
    %MinQValues=[];

    l=1;
    %Get optimal control policy as matrix
    for i=1:N1
        for j=1:N2
            for k=1:size(feasStates,3)
                for m=1:length(qValsMtx(i,j,k,:))
                    if (abs(min(qValsMtx(i,j,k,:))-qValsMtx(i,j,k,m))<epsilon4) && (feasStates(i,j,k)==1)
                        fullPolicyMtx(i,j,k,m)=1;
                    else
                        fullPolicyMtx(i,j,k,m)=0;
                    end
                end

                %ALSO, test if too MANY optimal optimal controls (more than one per state)...
                if(feasStates(i,j,k)==1)
                   %MinQValues=[MinQValues;i,j,k,min(qValsMtx(i,j,k,:))];
                   %NumOptCtrls_1=[NumOptCtrls_1;i,j,k,nnz(abs(qValsMtx(i,j,k,:)-min(qValsMtx(i,j,k,:)))<epsilon5)];
                   l=l+1;
                   NumOptCtrls_2=[NumOptCtrls_2;i,j,k,nnz(fullPolicyMtx(i,j,k,:))];
                end
            end
        end
    end
    
    numOptCtrls=sum(NumOptCtrls_2(:,4)); %Get TOTAL number of optimal ctrls
end