clearvars -except seqL MAX_LOAD MIN_LOAD E_MAX E_MIN D1Opt_State D2Opt_State

%Visualize all possible policies
E_Ind1=1:(E_MAX(1)-E_MIN(1)+1);
E_Ind2=1:(E_MAX(2)-E_MIN(2)+1);
indL=1:(MAX_LOAD-MIN_LOAD+1);

Dopt=D1Opt_State(:,:,:,1);
%Dopt=D2Opt_State(:,:,:,1);

figure
%%% Layered-Surfaces Visualization %%%%
for ind=1:length(indL)
    optD=Dopt(:,:,ind);
    Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
    surf(E_Ind2-1,E_Ind1-1,Z,optD)
    hold on
end
colorbar
xlabel('E2');ylabel('E1');zlabel('L');
title(['DOpt at t=9']);
hold off;

figure
%%% 3D Scatter Plot w/ Labels %%%
d=0.2; %displacement so the text does not overlay the data points
for ind=1:length(indL)
    [X,Y,Z]=Cuboid(E_Ind2-1, E_Ind1-1,ind-1);
    scatter3(X,Y,Z)
    numPts=length(E_Ind2)*length(E_Ind1)*1;
    optD=reshape(Dopt(:,:,ind),[numPts 1]);
    text(X,Y,Z+d,cellstr(num2str(optD)));
    hold on
end
xlabel('E2');ylabel('E1');zlabel('L');
title(['DOpt at t=9']);
hold off;