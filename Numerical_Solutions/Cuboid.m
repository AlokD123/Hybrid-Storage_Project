% - B.I SOLID POINT CLOUD CUBOID 
function [xp,yp,zp]=Cuboid(x,y,z) 
% - Generates coordinates for a solid cuboid composed of points

a=length(x);b=length(y);c=length(z); % assigning for coding simplicity 
% B.I.1. Finding all x coordinates 
xp=zeros(1,a*b*c); % preallocating 
for nx=1:a 
if nx==1,xp(1:c*b)=(repmat(x(nx),1,c*b)); 
else xp(((nx-1)*c*b)+1:nx*c*b)=repmat(x(nx),1,c*b);end 
end 
% B.I.2. Finding all y coordinates 
yp=zeros(1,b*c); % preallocating 
for ny=1:b 
if ny==1,yp(1:c)=repmat(y(ny),1,c); 
else yp(((ny-1)*c)+1:ny*c)=repmat(y(ny),1,c); 
end 
end 
yp=repmat(yp,1,a); 
% B.I.3. Finding all z coordinates 
zp=repmat(z,1,b*a); 
% B.I.4. Coordinates of the cuboid generated. 
%P=[xp;yp;zp]; 
end % To Optimize