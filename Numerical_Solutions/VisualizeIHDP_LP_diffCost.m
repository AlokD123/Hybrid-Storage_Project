%Script to visualize difference between (discretized) value functions

%PREREQUISITE: need to create matrix of difference between optimal values
%for subset of common states for OVFs with different sizes

INFCOST=0;

diffCosts=diffVal;

%diffCosts(isnan(diffCosts))=INFCOST;

%Padding (does not appear in 3D cost matrix)
diffCosts(end+1,:,:)=INFCOST;
diffCosts(:,end+1,:)=INFCOST;


    %Visualize all possible policies
    E_Ind1=1:(size(diffCosts,1)+1);
    E_Ind2=1:(size(diffCosts,2)+1);
    indL=1:(size(diffCosts,3)+1);

    %ti=0;tf=MAX_ITER-1;

    figure
    %%% Layered-Surfaces Visualization %%%%
    for t=1:1
        %subplot(1,1,t+1)
        for ind=1:size(diffCosts,3) %FULL LOAD NOT POSSIBLE, so only up to maximum possible
            costs=diffCosts(:,:,ind);
            Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
            surf(E_Ind2-1,E_Ind1-1,Z,costs)
            hold on
        end
        colorbar
        title(colorbar,'Difference')
        xlabel('Supercapacitor Energy (E2)');ylabel('Battery Energy (E1)');zlabel('Demand (L)');
        title('Difference in Optimal Value');
    end