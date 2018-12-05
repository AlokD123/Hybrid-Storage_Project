%Method to visualize optimal energy transfers
%v2: for realistic sizing of components

%% Parameter setup

E_SCALE_BATT=75/1000*2; %0.075 kWh/unit for battery
E_SCALE_SC_L=75;      %75 Wh/unit for supercapacitor and load
T_SCALE=2.5;          %2.5s/unit

%Create time axis
time=[1:length(optE1)]*T_SCALE;

%% Visualization

boolSingleSequence=1;
INFCOST=1e6;

clear Costs_LP;
Costs_LP=ConvCosts;

%Padding (does not appear in 3D cost matrix)
Costs_LP(end+1,:,:)=INFCOST;
Costs_LP(:,end+1,:)=INFCOST;


if(~boolSingleSequence)
    %Visualize all possible policies
    E_Ind1=1:(E_MAX(1)-E_MIN(1)+2);
    E_Ind2=1:(E_MAX(2)-E_MIN(2)+2);
    indL=1:(MAX_LOAD-MIN_LOAD+1);

    %ti=0;tf=MAX_ITER-1;

    figure
    %%% Layered-Surfaces Visualization %%%%
    for t=1:1
        %subplot(1,1,t+1)
        for ind=1:size(Costs_LP,3) %FULL LOAD NOT POSSIBLE, so only up to maximum possible
            costs=Costs_LP(:,:,ind);
            costs(costs>=INFCOST)=-100; %Replace all with infeasible costs for plotting
            Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
            surf(E_Ind2-1,E_Ind1-1,Z,costs)
            hold on
        end
        colorbar
        title(colorbar,'Optimal Value')
        xlabel('Supercapacitor Energy (E2)');ylabel('Battery Energy (E1)');zlabel('Demand (L)');
        title('Optimal Value Function for Discrete States');
    end

    
else
    %Visualize a possible sequence of loads, as previously found
    figure
    hold on;
    plot(time,optE1*E_SCALE_BATT,'-','MarkerSize',10);
    xlabel('Time (s)');
    xlim([0 max(time)]);
    ylabel('Battery energy (kWh)');
    ylim([0 inf])
    title('Variation in energy storage state under demand');
    yyaxis right
	plot(time,optE2*E_SCALE_SC_L,'-','MarkerSize',10); plot(time,[Load,0]*E_SCALE_SC_L,'--','MarkerSize',10);
    hold off;
    ylabel('Supercapacitor energy and demand (Wh)');
    ylim([-200 1000])
    legend('Battery (E1)','Supercapacitor (E2)','Demand (L)');
    
    %{
    figure
    hold on;
    plot(D1Opt,':*','MarkerSize',10); plot(D2Opt,':*','MarkerSize',10); plot(C2Opt,':*','MarkerSize',10); %; plot(Load,':o','MarkerSize',10);
    hold off;
    xlabel('Time (no units)');
    ylabel('Energy (no units)');
    title('Optimal policy vs. load');
    legend('Battery Discharge (D1)','Supercapacitor Discharge (D2)','Supercapacitor Charge (C2)'); %,'Demand (L)')
    axis([1 inf 0 max(MAX_DISCHARGE)]);
    %}
    
    %
    figure
    hold on;
    plot(U1Opt,':*','MarkerSize',10); plot(Load-U1Opt,':*','MarkerSize',10); plot(Load,':o','MarkerSize',10);
    hold off;
    xlabel('Time (no units)');
    ylabel('Energy (no units)');
    title('Optimal policy vs. load');
    legend('Battery Control (U1)','Supercapacitor Control (U2)','Demand (L)');
    axis([1 inf -max(MAX_CHARGE) max(MAX_DISCHARGE)]);
    %}
end