E_Ind1=1:(E_MAX(1)-E_MIN(1)+1);
E_Ind2=1:(E_MAX(2)-E_MIN(2)+1);
indL=1:(MAX_LOAD-MIN_LOAD+1);

nextE1=optNextE1(:,:,:,t);
nextE2=optNextE2(:,:,:,t);

figure
%%% Layered-Surfaces Visualization %%%%
for ind=1:length(indL)
    m=nextE1(:,:,ind);
    m(m==inf)=-10;
    Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
    surf(E_Ind2-1,E_Ind1-1,Z,m)
    hold on
end
colorbar
xlabel('E2');ylabel('E1');zlabel('L');
title(['Next State E1 at t=',num2str(t)]);
hold off;


figure
%%% Layered-Surfaces Visualization %%%%
for ind=1:length(indL)
    n=nextE2(:,:,ind);
    n(n==inf)=-10;
    Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
    surf(E_Ind2-1,E_Ind1-1,Z,n)
    hold on
end
colorbar
xlabel('E2');ylabel('E1');zlabel('L');
title(['Next State E2 at t=',num2str(t)]);
hold off;