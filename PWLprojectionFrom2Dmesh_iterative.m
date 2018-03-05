function [Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh_iterative(Xmn,Ymn,Omn,dt,nHalfSphereSamples)

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

% Create the angle vector in which the projections will be created
sphereAngs = linspace(0,pi,nHalfSphereSamples+1);
sphereAngs = sphereAngs(1:end-1);

% Piece wise linear sampling distance
ds = min([dx ; dy])/4;


% First create right hand and left hand side steps.
dtSampsRight = 0:dt:(ceil(R0/dt)+1)*dt;
dtSampsLeft = 0:-dt:-(ceil(R0/dt)+1)*dt;
% Concatenate excluding zeros
dtSamps = [dtSampsLeft(end:-1:2) dtSampsRight]; 

dsSamps = 0:ds:2*R0+ds;
% Iterate through projection angles, and create a corresponding set of 
% projections in that direction

% Create projection matrix_type
[T,S] = meshgrid(dtSamps,dsSamps);
pAngs = sphereAngs(:).';

Pij = zeros([numel(dtSamps) nHalfSphereSamples]);

for angIdx = 1:nHalfSphereSamples
  cXo = R0*cos(pAngs(angIdx)+pi/2);
  cYo = R0*sin(pAngs(angIdx)+pi/2);
  
  % Rotate
  cTr = bsxfun(@times,cos(pAngs(angIdx)+pi),T) + bsxfun(@times,-sin(pAngs(angIdx)+pi),S);
  cSr = bsxfun(@times,sin(pAngs(angIdx)+pi),T) + bsxfun(@times,cos(pAngs(angIdx)+pi),S);

  % Shift
  cTrs = bsxfun(@plus,cTr,cXo);
  cSrs = bsxfun(@plus,cSr,cYo);

%%%  %%%%%%%%%%%%%
%%%  % For debug %
%%%  %%%%%%%%%%%%%
%%%  plot(cTrs(1,:),cSrs(1,:),'.g','markersize',15);
%%%  hold on;
%%%  for j = 1:size(cTrs,2)
%%%    plot(cTrs(:,j),cSrs(:,j),'-b');
%%%  end
%%%  xlabel('x');
%%%  ylabel('y');
%%%  plot([minX minX maxX maxX minX],[minY maxY maxY minY minY],'-r','linewidth',2);
%%%  hold off;
%%%  pause(0.5);
%%%  end
%%%  %%%%%%%%%%%%%
%%%  % For debug %
%%%  %%%%%%%%%%%%%
  
  % Project
  cPij = interp2(Xmn,Ymn,Omn,cTrs(:),cSrs(:),'cubic',0);
  cPij = reshape(cPij,size(cTrs));

  % Perform trapezoidal integration
  cPij = 0.5*(cPij(2:end,:,:)+cPij(1:end-1,:,:))*ds;
  cPij = sum(cPij,1);

  % Permute for proper projection matrix
  Pij(:,angIdx) = cPij(:);
end



