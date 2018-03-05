function [Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh(Xmn,Ymn,Omn,dt,nHalfSphereSamples)
% [Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh(Xmn,Ymn,Omn,dt,nHalfSphereSamples)
% 
% Produces projections of a medium, while using a Piece-Wise Linear (PWL) 
% interpolation of the given medium samples. 
%
% Inputs:
% 
% Xmn -   [Ms,Ns] 'x' sample coordinates of 'Omn'.
% Ymn -   [Ms,Ns] 'y' sample coordinates of 'Omn'.
% Omn -   [Ms,Ns] Sampled values of the medium
% dt -    Sampling distance over 't'
% nHalfSphereSamples - Number of samples between phi = [0,pi)
%
% Outputs:
% Pij -   [M,nHalfSphereSamples] matrix denoting the 't' samples
%         along the rows and 'phi^i' samples along the columns.
% dtSamps - [1,M] Samples taken along 't' axis at 'dt' distances
% sphereAngs - [1,nHalfSphereSamples] samples taken along [0,pi)


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
ds = min([dx ; dy])/10;


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



% Perform 3rd dimension permuatation on the angles
pAngs = permute(sphereAngs(:),[3 2 1]);
Xo = R0*cos(pAngs+pi/2);
Yo = R0*sin(pAngs+pi/2);


% Rotate
Tr = bsxfun(@times,cos(pAngs+pi),T) + bsxfun(@times,-sin(pAngs+pi),S);
Sr = bsxfun(@times,sin(pAngs+pi),T) + bsxfun(@times,cos(pAngs+pi),S);

% Shift
Trs = bsxfun(@plus,Tr,Xo);
Srs = bsxfun(@plus,Sr,Yo);

%%%%%%%%%%%%%
% For debug %
%%%%%%%%%%%%%
%%for k = 1:size(Trs,3)
%%  plot(Trs(1,:,k),Srs(1,:,k),'.g','markersize',15);
%%  hold on;
%%  for j = 1:size(Trs,2)
%%    plot(Trs(:,j,k),Srs(:,j,k),'-b');
%%  end
%%  xlabel('x');
%%  ylabel('y');
%%  plot([minX minX maxX maxX minX],[minY maxY maxY minY minY],'-r','linewidth',2);
%%  hold off;
%%end
%%%%%%%%%%%%%
% For debug %
%%%%%%%%%%%%%

% Interpolate date for projection
Pij = interp2(Xmn,Ymn,Omn,Trs(:),Srs(:),'cubic',0);

Pij = reshape(Pij,size(Trs));
% Perform trapezoidal integration
Pij = 0.5*(Pij(2:end,:,:)+Pij(1:end-1,:,:))*ds;
Pij = sum(Pij,1);

% Permute for proper projection matrix
Pij = permute(Pij,[2 3 1]);
