function [optpath,totalcost  ] = dynamic_programming_shortestpath(cost,n  )
% Dynamic Programming for shortestpath problem
% This code has been written by W.Boonphakdee,D.Eng on Sept 20,2017
% If have any questions please contact me: warutboon@yahoo.com
% cost = matrix of network in form of an upper matrix
% cost=[ x x x 0 0 0
% 0 x 0 0 0 0
% 0 0 x x x 0
% 0 0 0 x x 0
% 0 0 0 0 x x
% 0 0 0 0 0 x ];
% n= number of node
% optpath= optimal path
% totalcost= cost of the optimal path
   %% Construct the functional equation matrix 
    fnmat=zeros(n,n);
    % Node n
    fn=zeros(n+1,1);
    for i=n:n
        for j=n:n
            fnmat(i,j)=cost(i,j);
            fn(i)=fnmat(i,j);
            iselect=i;
            jselect=j;
        end
    end
    % Node n-1 to 1
    for k=1:n-1
     i=n-k;
        for j=n:-1:i
            if cost(i,j)>0
               fnmat(i,j)=cost(i,j)+fn(j+1);
               iindex=i;
               jindex=j;
            end
              fn(iindex)=fnmat(iindex,jindex);
        end
         
        for j=n:-1:1
             if fnmat(i,j)>0 & fnmat(i,j)<=fn(i)
                iindex=i;
                jindex=j;
                fn(iindex)=fnmat(iindex,jindex);
             end
        end
    end
    
 %% Searching the  minimal functional equation matrix
minfn=zeros(n,n);
 for i=1:n
     for j=1:n
         
             if fn(i)==fnmat(i,j)
             imin=i;
             jmin=j;
             minfn(imin,jmin)=1;
             end
     end
 end
 %% Searching the optimal point matrix
 optpoint=zeros(n,n);
 % Node 1
 for i=1:1
     for j=1:n
         if minfn(i,j)==1
            iselect=i;
            jselect=j+1;
            optpoint(i,jselect)=1;
         end
     end
 end
 % Node 2 to node n
 for i=2:n
 
     for j= i:n
         if minfn(i,j)==1
             iselect=i;
             jselect=j+1;
             optpoint(iselect,jselect)=1;
         end
     end
 
 end
 %% Searching the optimal path
 
 
    % Find the number of optimal solution
   q=0;
   for j=1:n
       i=1;
       if optpoint(i,j)==1
           q=q+1;
       end;
   end;
 optpath=zeros(n,q);
 alltotalcost=sparse(q,1);
 for k=1:q
 
 for i=1:1
     for j=1:n+1
         if optpoint(i,j)==1
             iselect2=i;
             jselect2=j;
             optpath(jselect2,k)=cost(i,jselect2-1);
             optpoint(iselect2,jselect2)=0;
             break
         end
     end
 end
 
 for i=2:n
     for j=i:n+1
         if jselect2<n
            if optpoint(jselect2,j)==1 
               jselect3=j;
               optpath(jselect3,k)=cost(jselect2,jselect3-1);
              
               jselect2=jselect3;
            end
         end
     end
 end
 % Total cost
 
 totalcost=0;
 for i=1:n+1
     j=k;
     if optpath(i,j)>0
         totalcost=totalcost+optpath(i,j);
         optpath(n+2,k)=totalcost;
     end
    
 end
 alltotalcost(k,1)=totalcost;  
               
 end
             
end

