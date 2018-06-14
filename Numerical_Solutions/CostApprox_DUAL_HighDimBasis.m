%Script to TEST BEST APPROXIMATION, for one set of basis functions
%V3: Cost approximation for dual problem

n=4; %Order of polynomial approximation

%% Store FEASIBLE values of E1, E2 and L for each state

feasE1s=[]; feasE2s=[]; feasLs=[]; feasCtrls=[];

for p=1:P1*P2
    E_Ind_Vect_p=E_Ind_Vect{p};
    for i=1:length(E_Ind_Vect{p})
       if(i~=1 && E_Ind_Vect_p(i)==E_Ind_Vect_p(i-1))
           indL=indL+1;
       else
           indL=1;
       end
       E2=remainder(E_Ind_Vect_p(i),N2);
       E1=(E_Ind_Vect_p(i)-E2)/N2 + 1;
       L=indL;

       %Also add values to vectors, for each feasible state...
       feasE2s=[feasE2s;E2]; feasE1s=[feasE1s;E1]; feasLs=[feasLs;L]; %Add state
       feasCtrls=[feasCtrls;p]; %Add ctrl
    end
end


%% Test selected basis functions by 1-norm fitting
%Bases: all terms in order-n polynomial

Phi=[]; %Design matrix to create

%Adjoin feasible state-action vectors to form a ?x4 array
feasStatesArr=[feasE1s,feasE2s,feasLs,feasCtrls];

%Create design matrix with fitting functions up to order n
Phi=DesignMtx(feasStatesArr,d,n);

%Do 1-norm minimization on overdetermined system to find r...
cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable r_fit_dual(size(Phi,2))
    minimize( b'*Phi*r_fit_dual )
    subject to
        Q'*Phi*r_fit_dual == ones(size(Q',1),1)
        Phi*r_fit_dual >= 0
cvx_end

figure
hold on;
plot(Phi*r_fit_dual,Phi*r_fit_dual, '*');
plot(Phi*r_fit_dual,d, '*');
xlabel('Approximate Dual');
legend('Approximate Dual','Actual Dual');
title(strcat('Testing Dual Approximation Using Order-',num2str(n),' Fit (',num2str(length(r_fit_dual)),' Bases)'));