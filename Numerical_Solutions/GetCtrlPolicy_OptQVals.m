%Create matrix of optimal q-values and use to obtain optimal control policy
%qValsMtx(E1,E2,L,:) has q-values of state (E1,E2,L) for all controls

clear qValsMtx; global qValsMtx;
clear fullPolicyMtx; global fullPolicyMtx;

%N1=6; N2=5; M=10;
%P1=N1; P2=N2;

%Repeat for each decision (D1,D2 combination)
for p=1:p_max 
qSubVec=aug_Q((p-1)*N1*N2*(M-1)+1:p*N1*N2*(M-1));
%Format q-values sub-vector into E1xE2 matrices ((M-1) matrices, for each value of load)
    for i=0:(M-1)*N2:(M-1)*N2*(N1-1) 
        for j=0:(M-1)-1 
          for ind=(1+i+j):(M-1):((M-1)*(N2-1)+i+(j+1)) 
              if(mod(ind,(M-1)*N2)==0)
                ind2=((M-1)*N2-1-j)/(M-1)+1;
              else
                ind2=(mod(ind,(M-1)*N2)-1-j)/(M-1)+1;
              end
              qValsMtx(i/((M-1)*N2)+1,ind2,j+1,p)=qSubVec(ind);
          end
        end
    end
end

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
        end
    end
end