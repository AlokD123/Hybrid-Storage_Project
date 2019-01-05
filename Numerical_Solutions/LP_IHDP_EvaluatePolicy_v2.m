%GET INFINITE HORIZON POLICY USING LP SOLUTION
%V2: ADDED OPTIMAL POLICY INTERPOLATION

%NumIter=20; %Number of iterations of the policy to do

%E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
%E_MAX=[5;4]; %Maximum energy to be stored (upper bound)
%E1_INIT=E_MAX(1); 
%E2_INIT=E_MAX(2);

%STEP 1: initalize
optE1=[]; optE2=[]; D1Opt=[]; D2Opt=[]; Load=[];
global resL_Mult;

%STEP 2: evaluate policy for a random sequence of loads (ONLINE)
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads


NumIter=20; %Number of iterations of the policy to do

DeltaL_min=1/resL_Mult;

t_ind_VI=1; %Start evaluation

L=0; %Assume that the first demand is ZERO (starting smoothly)
while t_ind_VI<=NumIter
    %Set state index
    indE1=optE1(t_ind_VI)-E_MIN(1)+1;
    indE2=optE2(t_ind_VI)-E_MIN(2)+1;
    
    %% Create demand sequences and run online
    %Get optimal controls for given state
    [D1,D2]=GetPOpt(indE1,indE2,L);
    
    
    %While load makes for infeasible state, decrease till feasible
    while isempty(D1)||isempty(D2)||(D1==0&&D2==0&&L>0) %Infeasible if no optimal control OR 0 when non-zero load
        %SHOULD NOT GET HERE!!
        L=L-DeltaL_min;
        if(L<0)
            L=0;
            [D1,D2]=GetPOpt(indE1,indE2,L);
            break;
        end
        [D1,D2]=GetPOpt(indE1,indE2,L);
    end
    
    
    %Calculate the state these values of u and w will lead to, even if impossible...
    [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
    if abs(nextE2)<epsilon2 %Accounting for numerical errors
        nextE2=0;
    end
        
    %If leads to next state guaranteed out of bounds, decrease load
    while(nextE1>E_MAX(1)||nextE1<E_MIN(1)||nextE2>E_MAX(2)||nextE2<E_MIN(2)||(D1+D2-L)<-eps)
        L=L-DeltaL_min;%Decrement load
        %Increment out of bounds count
        countOOB=countOOB+1;
        if(L<0)
            L=0;
            %break;
        end
        
        %Calculate next state for new load value
        [D1,D2]=GetPOpt(indE1,indE2,L);
        [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    indL=(L-MIN_LOAD)*resL_Mult+1;
    Load(t_ind_VI)=L;          %Hold value of load (for reference)

    currL=L; %Store current load value
    
    %Get specific control sequence for given load sequence
    D1Opt(t_ind_VI)=D1;
    D2Opt(t_ind_VI)=D2;
    %Get optimal next state, starting at E_INIT
    optE1(t_ind_VI+1)=nextE1;
    optE2(t_ind_VI+1)=nextE2;
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        optE1(t_ind_VI+1)=optE1(t_ind_VI); optE2(t_ind_VI+1)=optE2(t_ind_VI);
        Load(t_ind_VI+1)=0;
    end
    
    %Finish iteration
    fbleCtrlOnGrd=[];
    fprintf("Completed iteration t=%d\n",t_ind_VI);
    
    t_ind_VI=t_ind_VI+1; %RUN ONLINE, next step
    
    
    %
    %% GENERATE NEXT DEMAND SAMPLE
    %**Note: GENERATE ACCORDING TO CUSTOM DISTRIBUTION (i.e. Binomial PDF)
    %Get number of feasible loads in next state
    numL_OffGrd=0; maxL_prev=0;
    nxtL=[];
    
    for D1_next=0:MAX_DISCHARGE(1)
        for D2_next=0:MAX_DISCHARGE(2)
            %Get number of FEASIBLE next loads
            %Check excess discharge condition
            if(~(D1_next>nextE1 || D2_next>nextE2))
                %For each perturbation at the NEXT time...
                for indL=(maxL_prev+1):((D1_next+D2_next-MIN_LOAD)*resL_Mult+1) %TRY ONLY PERTURBATIONS ABOVE PREVIOUS MAX
                    L=MIN_LOAD+(indL-1)/resL_Mult;
                    [next_nextE1,next_nextE2]=optNextStateLimited(nextE1,nextE2,D1_next,D2_next,L);
                    %Check other conditions
                    if(next_nextE1<=E_MAX(1) && next_nextE1>=E_MIN(1))
                        if(next_nextE2<=E_MAX(2) && next_nextE2>=E_MIN(2))
                            if(~((D1_next+D2_next-L)<0||(D1_next+D2_next-L)>MAX_CHARGE(2)))
                                %If feasible...
                                    nxtL=[nxtL;L]; %Store each load
                                    maxL_prev=maxL_prev+1;
                            end
                        end
                    end
                end
            end

        end
    end
    
    nxtL=unique(nxtL); %Remove duplicates
    numL_OffGrd=length(nxtL); %Get number of feasible loads
    
    %Create vector of next state probabilities for each of the next demands
    %Probabilities follow custom distribution sampled resL_Mult times more finely
    prob_nextL=ProbDistr_v2(numL_OffGrd,(currL-MIN_LOAD)*resL_Mult+1,(nxtL-MIN_LOAD)*resL_Mult+1,0);
    
    %SAMPLE FROM CONDITIONAL DISTRIBUTION TO GET INDEX OF NEXT LOAD in vector nxtL (NOT RELATIVE TO MIN_LOAD) 
    nxt_indL=find(mnrnd(1,prob_nextL),1);
    %Get demand value from index
    minL=nxtL(1);
    L=minL+(nxt_indL-1)/resL_Mult;
    %}
end