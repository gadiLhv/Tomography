clc
clear

close all

%%%%%%%%%%%%%%%%%%%%%%%
% Noise configuration %
%%%%%%%%%%%%%%%%%%%%%%%
SNR = 60;

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
sampDt = 0.1;
nHalfSphereSamples = 50;

nPts = [100,100];

Xrange = [-4,4];
Yrange = [-4,4];

circCenter = [-1.5,0.7];
circRad = 1;
circAmp = 0.4;


circCenter2 = [0.3,-1];
circRad2 = 2;
circAmp2 = 0.7;

squarePos = [0 2 ; 2 3];
squareAmp = 0.3;

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

[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh_iterative(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);
%[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);

% Determine variances to add SNR
Sigma_Pij = var(Pij);
SNR_lin = 10^(SNR/10);

Sigma_noise = Sigma_Pij/sqrt(SNR_lin);

Pij_noisy = Pij + bsxfun(@times,sqrt(Sigma_noise),randn(size(Pij)));


Omn_r = unfilteredBackprojection(Pij_noisy,sampDt,dtSamps,sphereAngs,Xmn,Ymn);
%Omn_r = unfilteredBackprojection_explicitFT(Pij,sampDt,dtSamps,sphereAngs,Xmn,Ymn);


figure('position',[70    200   1423    421]);
subplot(1,3,1);
imagesc(Xmn(1,:),Ymn(:,1),Omn);
xlabel('X');
ylabel('Y');
title('Original medium');
colorbar;
caxis([0 1]);
axis('square');

subplot(1,3,2);
imagesc(sphereAngs/pi,dtSamps,Pij);
xlabel('\phi/\Pi');
ylabel('t');
title('Projection data');
colorbar;
axis('square');

subplot(1,3,3);
imagesc(Xmn(1,:),Ymn(:,1),real(Omn_r));
xlabel('X');
ylabel('Y');
title('Reconstructed medium');
colorbar;
caxis([0 1]);
axis('square');

var(Omn_r(:)-Omn(:))

