%GET INFINITE HORIZON POLICY USING LP SOLUTION
%V4: COMBINED CONTROLS

%E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
%E_MAX=[5;4]; %Maximum energy to be stored (upper bound)
%E1_INIT=E_MAX(1); 
%E2_INIT=E_MAX(2);
global P0;

%Initialize
optE1=[]; optE2=[]; D1Opt=[]; D2Opt=[]; Load=[];

%%Evaluate policy for a random sequence of loads (ONLINE)
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads


NumIter=100000/100; %Number of iterations of the policy to do

DeltaL_min=1/resL_Mult;

t_ind_VI=1; %Start evaluation
L=0; %Assume that the first demand is ZERO (starting from rest)
while t_ind_VI<NumIter
    %Set state index
    indE1=RES_E1*(optE1(t_ind_VI)-E_MIN(1))+1;
    indE2=RES_E2*(optE2(t_ind_VI)-E_MIN(2))+1;
    
    %% Create demand sequences and run online
    %Get optimal controls for given state
    U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
    
    %While load makes for infeasible state, change till feasible. L=0
    %always feasible, so break if get there.
    while isempty(U1)||(L-U1)<-MAX_CHARGE(2)-eps||(L-U1)>MAX_DISCHARGE(2)+eps %Infeasible if no optimal control OR C2 constraints violated
        %SHOULD NOT GET HERE!!
        if L<0
            L=L+DeltaL_min;
            %L=0;                   %To test with no R.B.
        elseif L>0
            L=L-DeltaL_min;
            %L=max(0,L-DeltaL_min); %To test with no R.B.
        else
            U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
            break;
        end
        U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
    end
    
    
    %Calculate the state these values of u and w will lead to, even if impossible...
    [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
        
    %If leads to next state guaranteed out of bounds, change load
    while nextE1>E_MAX(1)+eps||nextE1<E_MIN(1)-eps||nextE2>E_MAX(2)+eps||nextE2<E_MIN(2)-eps||(L-U1)<-MAX_CHARGE(2)-eps||(L-U1)>MAX_DISCHARGE(2)+eps
        %Increment out of bounds count
        countOOB=countOOB+1;
        %Change load till feasible
        if L<0
            L=L+DeltaL_min;
            %L=0;                   %To test with no R.B.
        elseif L>0
            L=L-DeltaL_min;
            %L=max(0,L-DeltaL_min); %To test with no R.B.
        else            %For L=0, just get controls and exit
            U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
            [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
            break;
        end
        
        %Calculate next state for new load value
        U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
        [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    Load(t_ind_VI)=L;          %Hold value of load (for reference)

    currL=L; %Store current load value
    
    %Get specific control sequence for given load sequence
    U1Opt(t_ind_VI)=U1;
    %Get optimal next state, starting at E_INIT
    optE1(t_ind_VI+1)=nextE1;
    optE2(t_ind_VI+1)=nextE2;
    %Note: cross-check with optNextStateLtd?
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        optE1(t_ind_VI+1)=optE1(t_ind_VI); optE2(t_ind_VI+1)=optE2(t_ind_VI);
        Load(t_ind_VI+1)=0;
    end
    
    %Finish iteration
    fbleCtrlOnGrd=[];
    fprintf("Completed iteration t=%d\n",t_ind_VI);
    
    t_ind_VI=t_ind_VI+1; %RUN ONLINE, next step
    
    
    %% GENERATE NEXT DEMAND SAMPLE
    %**Note: GENERATE ACCORDING TO CUSTOM DISTRIBUTION (i.e. Binomial PDF)
    %Get number of feasible loads in next state
    numL_OffGrd=0; nxtL=[];
    
     for U1_Ind_next=1:(MAX_DISCHARGE(1)+MAX_CHARGE(1))*RES_U1+1
        U1_next=(U1_Ind_next-1)/RES_U1-MAX_CHARGE(1);
        %Get number of FEASIBLE next loads...
        
        %Check excess discharge condition
        if( ~(U1_next>nextE1 || U1_next<(nextE1-E_MAX(1))) )
            %For each perturbation at the NEXT time...
            for L=linspace(-MAX_CHARGE(2)+U1_next,MAX_DISCHARGE(2)+U1_next,resL_Mult*(MAX_DISCHARGE(2)+MAX_CHARGE(2))+1)
                [next_nextE1,next_nextE2]=optNextStateLimited_v3(nextE1,nextE2,U1_next,L);
                %Check other conditions
                if(next_nextE1<=E_MAX(1) && next_nextE1>=E_MIN(1))
                    if(next_nextE2<=E_MAX(2) && next_nextE2>=E_MIN(2))
                        if(~( (-MAX_CHARGE(2)+U1_next)>L || L>(MAX_DISCHARGE(2)+U1_next) ))
                            %If feasible...
                            nxtL=[nxtL;L]; %Store each load
                        end
                    end
                end
            end
        end
     end
    
    %Get indices 
    idxCurrL=round((currL-MIN_LOAD)*resL_Mult+1);
    idxNxtL=round((nxtL-MIN_LOAD)*resL_Mult+1);
     
    idxNxtL=unique(idxNxtL); %Remove duplicates
    numL_OffGrd=length(idxNxtL); %Get number of feasible loads
    
    %Create vector of next state probabilities for each of the next demands
    %Probabilities follow custom distribution sampled resL_Mult times more finely
    prob_nextL=ProbDistr_v2(numL_OffGrd,idxCurrL,idxNxtL,0);
    
    %SAMPLE FROM CONDITIONAL DISTRIBUTION TO GET INDEX OF NEXT LOAD in vector nxtL (NOT RELATIVE TO MIN_LOAD) 
    nxt_indL=find(mnrnd(1,prob_nextL),1);
    %Get demand value from index
    minL=nxtL(1);
    L=minL+(nxt_indL-1)/resL_Mult;
end