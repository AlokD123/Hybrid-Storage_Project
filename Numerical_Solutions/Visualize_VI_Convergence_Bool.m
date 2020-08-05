%Script to verify iterative minimization of value function (with boolean indicators of decrease in each state)
%clearvars -except seqL MAX_LOAD MIN_LOAD E_MAX E_MIN BOOL_VI_CONV BOOL_VI_CONV_PREV

E_Ind1=1:(E_MAX(1)-E_MIN(1)+1);
E_Ind2=1:(E_MAX(2)-E_MIN(2)+1);
indL=1:(MAX_LOAD-MIN_LOAD+1);

bool=BOOL_VI_CONV(:,:,:);

figure
%%% Layered-Surfaces Visualization %%%%
for ind=1:length(indL)
    b=bool(:,:,ind);
    Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
    surf(E_Ind2-1,E_Ind1-1,Z,b)
    hold on
end
colorbar
xlabel('E2');ylabel('E1');zlabel('L');
title(['Convergence at t=',num2str(t)]);
hold off;


figure
%%% Layered-Surfaces Visualization %%%%
for ind=1:length(indL)
    cost=V(:,:,ind,t);
    cost(cost==inf)=-10;
    Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
    surf(E_Ind2-1,E_Ind1-1,Z,cost)
    hold on
end
colorbar
xlabel('E2');ylabel('E1');zlabel('L');
title(['Cost at t=',num2str(t)]);
hold off;

% figure
% %%% 3D Scatter Plot w/ Labels %%%
% d=0.2; %displacement so the text does not overlay the data points
% for ind=1:length(indL)
%     [X,Y,Z]=Cuboid(E_Ind2-1, E_Ind1-1,ind-1);
%     scatter3(X,Y,Z)
%     numPts=length(E_Ind2)*length(E_Ind1)*1;
%     b=reshape(bool(:,:,ind),[numPts 1]);
%     text(X,Y,Z+d,cellstr(num2str(b)));
%     hold on
% end
% xlabel('E2');ylabel('E1');zlabel('L');
% title(['Convergence at t=',num2str(t)]);
% hold off;