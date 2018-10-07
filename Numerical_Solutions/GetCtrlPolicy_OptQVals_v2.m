%Create matrix of optimal q-values and use to obtain optimal control policy
%qValsMtx(E1,E2,L,:) has q-values of state (E1,E2,L) for all controls

global qValsMtx; global fullPolicyMtx;
qValsMtx=[]; fullPolicyMtx=[];

%N1=6; N2=5; M=10;
%P1=N1; P2=N2;

%Define Sigma as the total number of possible loads over all E-states and
%controls
Sigma=size(E_Ind_MtxALL,2); %i.e. -E_MAX -> M-2


%Create CONCATENATED matrix of feasible controls per state
feasStates_ctrl_i=[];
feasStatesArr=[];
feasStatesCatArr=[];
feasStates_ctrl={};

for i=1:length(feasStatesArr_ctrl)
    feasStatesArr=[feasStatesArr;feasStatesArr_ctrl{i}];
    for j=1:size(feasStatesArr_ctrl{i},1)
        feasStates_ctrl_i(feasStatesArr_ctrl{i}(j,1),feasStatesArr_ctrl{i}(j,2),feasStatesArr_ctrl{i}(j,3))=1;
    end
    feasStates_ctrl{1,i}=padarray(feasStates_ctrl_i,[N1-size(feasStates_ctrl_i,1),N2-size(feasStates_ctrl_i,2),Sigma-size(feasStates_ctrl_i,3)],'post');
    feasStates_ctrl_i=[];
    feasStatesCatArr(:,:,:,i)=feasStates_ctrl{i}; %Same size as qValsMtx
end


%Create Q-value matrix
for p=1:(p_max-CountInfCtrls) 
    qSubVec=aug_Q((p-1)*N1*N2*Sigma+1:p*N1*N2*Sigma);
    %Format q-values sub-vector into E1xE2 matrices (Sigma matrices, one for each value of load)
    for i=0:Sigma*N2:Sigma*N2*(N1-1) 
        for j=0:Sigma-1 
          for ind=(1+i+j):Sigma:(Sigma*(N2-1)+i+(j+1)) 
              if(mod(ind,Sigma*N2)==0)
                ind2=(Sigma*N2-1-j)/Sigma+1;
              else
                ind2=(mod(ind,Sigma*N2)-1-j)/Sigma+1;
              end
              qValsMtx(i/(Sigma*N2)+1,ind2,j+1,p)=qSubVec(ind);
          end
        end
    end
end

qValsMtx(qValsMtx==0)=Inf;


%Get optimal control policy as matrix
for i=1:N1
    for j=1:N2
        for k=1:Sigma
            for m=1:length(qValsMtx(i,j,k,:))
                if (abs(min(qValsMtx(i,j,k,:))-qValsMtx(i,j,k,m))<epsilon4) && feasStatesCatArr(i,j,k,m)==1 % ismember([i+E_MIN(1)-1,j+E_MIN(2)-1,k+MIN_LOAD-1],feasStatesArr_ctrl{m},'rows')
                    fullPolicyMtx(i,j,k,m)=1;
                else
                    fullPolicyMtx(i,j,k,m)=0;
                end
            end
        end
    end
end

%Find any states that do not have optimal controls
missingStatesMtx=~(max(fullPolicyMtx,[],4)+~feasStates);

%For each of these states, set optimal control to be that giving NEXT
%LOWEST Q-value amongst controls in feasible set
for i=1:N1
    for j=1:N2
        for k=1:Sigma
            if missingStatesMtx(i,j,k)==1
                feasOptCtrlIdxs=find(feasStatesCatArr(i,j,k,:)==1);
                optCtrlIdx=feasOptCtrlIdxs(1);
                tempMinQ=qValsMtx(i,j,k,feasOptCtrlIdxs(1));
                for m=2:length(feasOptCtrlIdxs)
                    if qValsMtx(i,j,k,feasOptCtrlIdxs(m))<tempMinQ
                        tempMinQ=qValsMtx(i,j,k,feasOptCtrlIdxs(m));
                        optCtrlIdx=feasOptCtrlIdxs(m);
                    end
                end
                fullPolicyMtx(i,j,k,optCtrlIdx)=1; %Set this to be optimal control
            end
        end
    end
end