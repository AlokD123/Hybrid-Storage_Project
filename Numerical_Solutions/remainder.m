function [ E_Ind2 ] = remainder( E_Ind, N2 )
%remainder: returns index associated with state E2 for given linear indexing of
%both states

if (rem(E_Ind,N2)==0)
   E_Ind2=N2;
else
   E_Ind2=rem(E_Ind,N2);
end

end

