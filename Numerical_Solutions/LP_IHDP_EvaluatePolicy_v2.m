%GET INFINITE HORIZON POLICY USING LP SOLUTION
%V2: ADDED OPTIMAL POLICY INTERPOLATION

%NumIter=20; %Number of iterations of the policy to do

%E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
%E_MAX=[5;4]; %Maximum energy to be stored (upper bound)
%E1_INIT=E_MAX(1); 
%E2_INIT=E_MAX(2);

%STEP 1: initalize
optE1=[]; optE2=[]; D1Opt=[]; D2Opt=[]; Load=[];
resL_Mult=2; %Set resolution of demands by a dividing factor (natural number)

%STEP 2: evaluate policy for a random sequence of loads (ONLINE)
%Set up matricee
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads

%OLD - load sequence
%seqL=[]; %Reset load sequence
%seqL=ones(NumIter,1); %CONSTANT LOAD
%seqL=[3 1 0 1 0 2 0 1 1 0 3 1 0 1 0 2 0 1 1 0]; %FLUCTUATING LOAD
%seqL=[1 1 0 0 1 1 2 2 1 0 1 1 0 0 1 1 2 2 1 0]; %SMOOTH LOAD

%FLUCTUATING LOAD
%{
seqL=[0.25 1 0.75 0 0.5 1 1.75 2 1.25 2 2.5 2.75 1.5 3 2.25 2.75 3.5 3 3.25 3.25]; %FLUCTUATING LOAD
for i=1:1
    seqL=[seqL,2.75*(1+normrnd(0.2,0.1,1,140)).*(1+0.1*sin(linspace(1,50,140))).*(1-0.2*sin(linspace(1,10,140)))];

end
%}

%SMOOTH LOAD
%{
seqL=[0.25 0 0.25 0.5 0.25 0.75 1.25 0.5 1.75 2 2.25 1.75 2.25 2.75 2.75 3 2.75 3 3.25 3.5];
for i=1:1
    seqL=[seqL,fliplr(seqL(end/2:end)).*linspace(1,1.75,length(seqL(end/2:end)))];
    seqL=[seqL,(1+normrnd(0.2,0.1,1,129)).*linspace(3.2,1,129).*(1+0.1*sin(linspace(1,10,129)))];
end
%}

%RAMP LOAD
%{
seqL=[0.25 0 0.25 0.5 0.75 1 1 0.75 1 1.25 1.25 1.25 1 1.5 1.5 1.75 2 2.25 2.25 2.25];
for i=1:1
    seqL=[seqL,2*rand(1)*seqL(end/2:end)+0.5];
    seqL=[seqL,rand(1)*fliplr(seqL(end/2:end))];
    seqL=[seqL,rand(1)*seqL(end/2:end)/2];
    seqL=[seqL,rand(1)*fliplr(seqL(1:end))/2];
    seqL=[seqL,fliplr(seqL(end-17:end))+0.25];
end
%}

%seqL=ones(1,10);

%{
seqL=seqL*3/4;
seqL(seqL>3)=3;
seqL(seqL<0)=0;


plot(seqL);
ylim([0 8])
%}

NumIter=length(seqL); %Number of iterations of the policy to do

DeltaL_min=0.2;

t_ind_VI=1; %Start evaluation

L=0; %Assume that the first demand is ZERO (starting smoothly)
while t_ind_VI<=NumIter
    %Set state index
    indE1=optE1(t_ind_VI)-E_MIN(1)+1;
    indE2=optE2(t_ind_VI)-E_MIN(2)+1;
    
    %% OLD- create demand sequence
    
    %{
    %OPTION 1: UNUSED
    %Create random demand from IID Uniform probability sequence
    MAX_LOAD_STATE=optE1(t_ind_VI)+optE2(t_ind_VI)-1; %Maximum possible load limited to total energy stored in that state
    if(MAX_LOAD_STATE==Inf)
        MAX_LOAD_STATE=MAX_LOAD;
    end
    %randL=randi(MAX_LOAD_STATE-MIN_LOAD+1,1,1)+MIN_LOAD-1;

    %Option 2: create sample sequence of pseudo-random demands
    %Select between sequences
    %L=randL;
    L= seqL(t_ind_VI);
    %L=min(t_ind_VI,MAX_LOAD_STATE); %Ramp
    
    %If demand is infeasible for current E-state, reduce HERE
    %When generating truly randomly using transition probabilities, this is automatically enforced
    %1) Case: load is too high for E-state
    if L>(E_MAX(1)+E_MAX(2))
        L=max(0,(E_MAX(1)+E_MAX(2))-( L-(E_MAX(1)+E_MAX(2) )) );
    end
    %2) Case: no control in feasible set leads to a next E-state that is feasible
    boolFeasL=0;
    %Decrease load till feasible for some set of controls on grid
    while L>=0 && boolFeasL==0
        %Check if feasible for these controls
        for D1=0:MAX_DISCHARGE(1) 
            for D2=0:MAX_DISCHARGE(2)
                if D1+D2>=L %For each control on grid that is feasible...
                    [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);

                    %If the next E-state is feasible, this is a feasible load (on or OFF grid)
                    if ~(nextE1>E_MAX(1)||nextE1<E_MIN(1)||nextE2>E_MAX(2)||nextE2<E_MIN(2))
                        %Store (one) feasible control on grid, to confirm ones off grid are unnecessarily infeasible
                        fbleCtrlOnGrd=[D1,D2,L];
                        %Exit loop with this value of L (feasible)
                        D1=MAX_DISCHARGE(1)+1; D2=MAX_DISCHARGE(2)+1;
                        boolFeasL=1;
                    end
                end
            end
        end
        if boolFeasL==0 %If infeasible, decrease load
           L=L-DeltaL_min; 
        end
    end
    
    %If supercap is empty and load is more than D1-max, reduce to range in [0, D1-max]
    if abs(indE2-1)<epsilon2 && L>MAX_DISCHARGE(1)
        L=max(0,MAX_DISCHARGE(1)-(L-MAX_DISCHARGE(1)));
    end
    
    %}
    
    %% Run sequences online
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
    indL=L-MIN_LOAD+1;
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
                for indL=(maxL_prev+1):(D1_next+D2_next-MIN_LOAD+1) %TRY ONLY PERTURBATIONS ABOVE PREVIOUS MAX
                    L=indL+MIN_LOAD-1;
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
    %First, assume uniform
    %prob_nextL=1/numL_OffGrd*ones(1,numL_OffGrd)';
    %Then, transform to custom distribution sampled resL_Mult times more finely
    prob_nextL=ProbDistr_v2(numL_OffGrd,currL-MIN_LOAD+1,nxtL-MIN_LOAD+1,resL_Mult);
    
    %SAMPLE FROM CONDITIONAL DISTRIBUTION TO GET INDEX OF NEXT LOAD in vector nxtL (NOT RELATIVE TO MIN_LOAD) 
    nxt_indL=find(mnrnd(1,prob_nextL),1);
    %Get demand value from index
    minL=nxtL(1);
    L=minL+nxt_indL/resL_Mult;
    %}
end