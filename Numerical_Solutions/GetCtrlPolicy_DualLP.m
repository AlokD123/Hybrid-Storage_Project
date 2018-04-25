clear policyMtx;

N1=6; N2=5; M=10;
P1=N1; P2=N2;

%Repeat for each decision (D1,D2 combination)
for(p=1:P1*P2)
piSubVec=pi((p-1)*N1*N2*M+1:p*N1*N2*M);
%Format binary policy sub-vector (for cases when applying the single decision) into E1xE2 matrices (M matrices, for each value of load)
    for(i1=0:M*N2:M*N2*(N1-1))
        for(j=0:M-1)
          for(ind=(1+i1+j):M:(M*(N2-1)+i1+(j+1)))
              if(mod(ind,M*N2)==0)
                ind2=(M*N2-1-j)/M+1;
              else
                ind2=(mod(ind,M*N2)-1-j)/M+1;
              end
              policyMtx(i1/(M*N2)+1,ind2,j+1,p)=piSubVec(ind);
          end
        end
    end
end