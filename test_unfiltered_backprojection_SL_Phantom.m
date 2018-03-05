clc
clear

close all

%%%%%%%%%%%%%%%%%%
%% Debug config %%
%%%%%%%%%%%%%%%%%%
%
%% Sampling distance along line. 
%sampDt = 0.5;
%nHalfSphereSamples = 10;
%
%
% Medium mesh config
nPts = [100,100];
Xrange = [-4,4];
Yrange = [-4,4];

%circCenter = [-1.5,0.7];
%circRad = 1;
%circAmp = 1.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dummy data parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sampling distance along line. 
sampDt = 0.1;
nHalfSphereSamples = 50;


Xrange = [-4,4];
Yrange = [-4,4];

Omn = phantom('Modified Shepp-Logan',100);

nPts = size(Omn);

% Create mesh
[Xmn,Ymn] = meshgrid(linspace(Xrange(1),Xrange(2),nPts(1)),linspace(Yrange(1),Yrange(2),nPts(2)));

if(exist('show_slices.mat') == 0)

%% Create set of projections

[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh_iterative(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);


%Omn_r = unfilteredBackprojection(Pij,sampDt,dtSamps,sphereAngs,Xmn,Ymn);
%Omn_r = unfilteredBackprojection_explicitFT(Pij,sampDt,dtSamps,sphereAngs,Xmn,Ymn);
Omn_r = unfilteredBackprojection_TD(Pij,sampDt,dtSamps,sphereAngs,Xmn,Ymn);

end

load show_slices.mat
system('rm show_slices.mat');

figure('position',[70    200   1520    520]);
subplot(1,3,1);
origHdl = imagesc(Xmn(1,:),Ymn(:,1),Omn);
xlabel('X');
ylabel('Y');
title('Original medium');
colorbar;
caxis([0 1]);
axis('square');

subplot(1,3,2);
Pij_curr = NaN*ones(size(Pij));
slicesHdl = imagesc(sphereAngs/pi,dtSamps,Pij_curr);
xlabel('\phi/\Pi');
ylabel('t');
title('Projection data');
colorbar;
axis('square');

subplot(1,3,3);
cRec = zeros(size(Omn));
recHdl = imagesc(Xmn(1,:),Ymn(:,1),real(cRec));
xlabel('X');
ylabel('Y');
title('Reconstructed medium');
colorbar;
caxis([0 1]);
axis('square');

for sliceIdx = 1:size(layered_Omn,2)
  Pij_curr(:,sliceIdx) = Pij(:,sliceIdx);
  cRec = cRec + reshape(layered_Omn(:,sliceIdx),size(Omn));
  set(slicesHdl,'CData',Pij_curr);
  set(recHdl,'CData',real(cRec));
  
  pause(3/nHalfSphereSamples);
end

