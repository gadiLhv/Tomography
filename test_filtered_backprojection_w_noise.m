clc
clear

close all

%%%%%%%%%%%%%%%%%%%%%%%
% Noise configuration %
%%%%%%%%%%%%%%%%%%%%%%%
SNR = 30;
filterBW = 0.25;   % Sets the -3dB point of the window in the 20% of the available bandwidth

%%%%%%%%%%%%%%%%%%
%% Debug config %%
%%%%%%%%%%%%%%%%%%
%
%% Sampling distance along line. 
%sampDt = 0.5;
%nHalfSphereSamples = 10;
%
%
%% Medium mesh config
%nPts = [10,10];
%Xrange = [-4,4];
%Yrange = [-4,4];
%
%circCenter = [-1.5,0.7];
%circRad = 1;
%circAmp = 1.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dummy data parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sampling distance along line. 
sampDt = 0.025;
nHalfSphereSamples = 200;

nPts = [100,100];

Xrange = [-4,4];
Yrange = [-4,4];

circCenter = [-1.5,0.7];
circRad = 1;
circAmp = 0.1;

circCenter2 = [0.3,-1];
circRad2 = 2;
circAmp2 = 0.05;

squarePos = [0 2 ; 2 3];
squareAmp = 0.075;

rimWidth = 3;
rimAmp = 2;

% Create mesh
[Xmn,Ymn] = meshgrid(linspace(Xrange(1),Xrange(2),nPts(1)),linspace(Yrange(1),Yrange(2),nPts(2)));

% Create the amplitude data
Omn = zeros(size(Xmn));
binCirc = sqrt((Xmn-circCenter(1)).^2+(Ymn-circCenter(2)).^2) <= circRad;
%Omn(binCirc) = circAmp*(circRad-sqrt((Xmn(binCirc)-circCenter(1)).^2+(Ymn(binCirc)-circCenter(2)).^2))/circRad;
Omn(binCirc) = circAmp;

binCirc2 = sqrt((Xmn-circCenter2(1)).^2+(Ymn-circCenter2(2)).^2) <= circRad2;
Omn(binCirc2) = circAmp2*(sqrt((Xmn(binCirc2)-circCenter2(1)).^2+(Ymn(binCirc2)-circCenter2(2)).^2))/circRad2;;

Omn(Xmn >= squarePos(1,1) & Xmn <= squarePos(2,1) & Ymn >= squarePos(1,2) & Ymn <= squarePos(2,2)) = squareAmp;

%binMap = true(size(Omn));
%binMap(rimWidth:end-rimWidth,rimWidth:end-rimWidth) = ~binMap(rimWidth:end-rimWidth,rimWidth:end-rimWidth);
%Omn(binMap) = rimAmp;

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


