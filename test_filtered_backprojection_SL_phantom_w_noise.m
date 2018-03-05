clc
clear

close all

%%%%%%%%%%%%%%%%%%%%%%%
% Noise configuration %
%%%%%%%%%%%%%%%%%%%%%%%
SNR = 30;
filterBW = 0.25;   % Sets the -3dB point of the window in the 20% of the available bandwidth


% Sampling distance along line. 
sampDt = 0.025;
nHalfSphereSamples = 200;

Xrange = [-4,4];
Yrange = [-4,4];

Omn = phantom(100);

nPts = size(Omn);

% Create mesh
[Xmn,Ymn] = meshgrid(linspace(Xrange(1),Xrange(2),nPts(1)),linspace(Yrange(1),Yrange(2),nPts(2)));


%% Create set of projections

%[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);
[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh_iterative(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);

% Determine variances to add SNR
Sigma_Pij = var(Pij);

SNR_lin = 10^(SNR/10);

Sigma_noise = Sigma_Pij/sqrt(SNR_lin);

% Noise
Nij = bsxfun(@times,sqrt(Sigma_noise),randn(size(Pij)));
% Add noise
Pij_noisy = Pij + Nij;


Omn_r = filteredBackprojection(Pij_noisy,sampDt,filterBW,dtSamps,sphereAngs,Xmn,Ymn);


figure('position',[70    200   1423    421]);
subplot(1,3,1);
imagesc(Xmn(1,:),Ymn(:,1),Omn);
xlabel('X','fontsize',17);
ylabel('Y','fontsize',17);
title('Original medium','fontsize',20);
hdl = colorbar;
set(hdl,'fontsize',17);
caxis([0 max(Omn(:))]);
axis('square');

set(gca,'fontsize',17);

subplot(1,3,2);
imagesc(sphereAngs/pi,dtSamps,Pij);
xlabel('\phi/\Pi','fontsize',17);
ylabel('t','fontsize',17);
title('Projection data','fontsize',20);
hdl = colorbar;
set(hdl,'fontsize',17);
axis('square');

set(gca,'fontsize',17);

subplot(1,3,3);
imagesc(Xmn(1,:),Ymn(:,1),real(Omn_r));
xlabel('X','fontsize',17);
ylabel('Y','fontsize',17);
title('Reconstructed medium','fontsize',20);
hdl = colorbar;
set(hdl,'fontsize',17);
caxis([0 max(Omn(:))]);
axis('square');

set(gca,'fontsize',17);


