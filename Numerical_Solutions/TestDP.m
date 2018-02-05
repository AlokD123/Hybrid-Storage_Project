%%TESTING: evaluate policy for a random sequence of loads
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
%Find total initial storage (for creating sample load sequences)
TOT_INIT_E=optE1(1)+optE2(1);
% Count errors, for DEBUGGING
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads
while t<=(LAST_ITER-1)
    %Set state index
    indE1=optE1(t)-E_MIN(1)+1;
    indE2=optE2(t)-E_MIN(2)+1;
    %Create random demand in IID Uniform probability sequence
    MAX_LOAD_STATE=optE1(t)+optE2(t)-1; %Maximum possible load limited to total energy stored in that state
    randL=randi(MAX_LOAD_STATE-MIN_LOAD+1)+MIN_LOAD-1;
    %Or, use sample sequence of pseudo-random demands
    %seqL=max(randi(2)-1,1); %TOT_INIT_E-(t+2)+
    %Select between sequences
    %L=randL;
    L=max(TOT_INIT_E-(t+2),0);
    %L=seqL(t);
    %If leads to next state guaranteed out of bounds, decrease load
    while(optNextE1(indE1,indE2,L-MIN_LOAD+1,t)==inf || optNextE2(indE1,indE2,L-MIN_LOAD+1,t)==inf)
        L=L-1;
        countOOB=countOOB+1;
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    indL=L-MIN_LOAD+1;
    Load(t)=L;          %Hold value of load (for reference)
    %Get optimal costs for this sequence of loads
    optV(t)=V(indE1,indE2,indL,t);
    %Get costs of other possible states for this state and load??
    %(For comparison)

    %Get optimal policy for this sequence of loads
    D1Opt(t)=D1Opt_State(indE1,indE2,indL,t);
    D2Opt(t)=D2Opt_State(indE1,indE2,indL,t);
    %Get optimal next state, starting at E_INIT
    optE1(t+1)=optNextE1(indE1,indE2,indL,t);
    optE2(t+1)=optNextE2(indE1,indE2,indL,t);
    %Note: cross-check with optNextStateLtd?
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        D1Opt(t+1)=0; D2Opt(t+1)=0;
        optE1(t+1)=optE1(t); optE2(t+1)=optE2(t);
        optV(t+1)=V(optE1(t+1)-E_MIN(1)+1,optE2(t+1)-E_MIN(2)+1,indL,t);
        Load(t+1)=0;
        t=LAST_ITER;
    end
    t=t+1;
end