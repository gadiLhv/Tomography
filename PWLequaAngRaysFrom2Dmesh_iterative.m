function [Rij,dgSamps,D,b] = PWLequaAngRaysFrom2Dmesh_iterative(Xmn,Ymn,Omn,dg,gMax,D,nSphereSamples)

% Define object domain sphere
minX = min(Xmn(:));
maxX = max(Xmn(:));
minY = min(Ymn(:));
maxY = max(Ymn(:));

% Determine what is the squre size of the data mesh. Sampled points are 
% considered the square center.
dx = Xmn(1,2)-Xmn(1,1);
dy = Ymn(2)-Ymn(1);


% Determine object domain center
X0 = mean(Xmn(1,:).');
Y0 = mean(Ymn(:,1));

% Set object domain radius (bounding cirle radius)
R0 = sqrt((maxX-minX).^2+(maxY-minY).^2)/2;
D = max([R0 ; D]);


% Create the angle vector in which the projections will be created
b = linspace(0,2*pi,nSphereSamples+1);
b = b(1:end-1);

% Piece wise linear sampling distance
ds = min([dx ; dy])/4;


% First create right hand and left hand side steps.
dgSampsRight = 0:dg:(ceil(gMax/dg)+1)*dg;
dgSampsLeft = 0:-dg:-(ceil(gMax/dg)+1)*dg;
% Concatenate excluding zeros
dgSamps = [dgSampsLeft(end:-1:2) dgSampsRight]; 

% Samples along the ray emanating from each origin

dsSamps = 0:ds:(R0+D+ds);
% Iterate through projection angles, and create a corresponding set of 
% projections in that direction

% Create projection matrix_type
[G,S] = meshgrid(dgSamps,dsSamps);

Rij = zeros([size(G,2) nSphereSamples]);

for angIdx = 1:nSphereSamples
  cXo = D*cos(b(angIdx)+pi/2);
  cYo = D*sin(b(angIdx)+pi/2);
  
  
  % Rotate first to align with local projection point coordinate
  % system, in which the pointer towards the global system
  % is the "y" axis (or "s")
  cXr = cos(pi/2-G).*S;
  cYr = sin(pi/2-G).*S;

  % Now shift it so the local coordinate system will be aligned
  % with the global one.
  cXrs = cXr;
  cYrs = cYr - D;
  
  % Now, rotate again to go back to global coordinate system
  cXrsr =  cXrs.*cos(pi+b(angIdx)) - cYrs.*sin(pi+b(angIdx));
  cYrsr =  cXrs.*sin(pi+b(angIdx)) + cYrs.*cos(pi+b(angIdx));
  
%%  % plot (for debug)
%%%  plot(cXo,cYo,'.r','markersize',25);
%%  plot([minX minX maxX maxX minX],[minY maxY maxY minY minY],'-g','linewidth',5);
%%  hold on;
%%  axis([-(D+R0) (D+R0) -(D+R0) (D+R0)]*1.01);
%%  for idx = 1:size(cXrsr,2)
%%    plot(cXrsr(:,idx),cYrsr(:,idx),'-b');
%%    plot(cXrsr(1,idx),cYrsr(1,idx),'.g','markersize',20);
%%    plot(cXrsr(end,idx),cYrsr(end,idx),'.r','markersize',20);
%%    pause(0.01);
%%  end
%%  hold off;
%%  axis('square');
%%  pause(0.5);

  % Project
  cRij = interp2(Xmn,Ymn,Omn,cXrsr(:),cYrsr(:),'cubic',0);
  cRij = reshape(cRij,size(cXrsr));
  
  % Perform trapezoidal integration
  cRij = 0.5*(cRij(2:end,:,:)+cRij(1:end-1,:,:))*ds;
  cRij = sum(cRij,1);
  
  % Permute for proper projection matrix
  Rij(:,angIdx) = cRij(:);
end

dgSamps = dgSamps(:);



