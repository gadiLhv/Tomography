clc
clear

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dummy data parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sampling distance along line. 
sampDt = 0.1;
nHalfSphereSamples = 10;

nPts = 100;

Xrange = [-4,4];
Yrange = [-4,4];


% Create mesh
[Xmn,Ymn] = meshgrid(linspace(Xrange(1),Xrange(2),nPts),linspace(Yrange(1),Yrange(2),nPts));

Omn = phantom(nPts);

%% Create set of projections

[Pij,dtSamps,sphereAngs] = PWLprojectionFrom2Dmesh(Xmn,Ymn,Omn,sampDt,nHalfSphereSamples);

% Reconstruct
Omn_r = unfilteredBackprojection(Pij,sampDt,dtSamps,sphereAngs,Xmn,Ymn);

% Display results
fHdl = figure('position',[70    200   1423    421]);
subplot(1,3,1);
imagesc(Xmn(1,:),Ymn(:,1),Omn);
xlabel('\bf{x}','fontsize',14);
ylabel('\bf{y}','fontsize',14);
title('\bf{Original medium}','fontsize',14);
hdl = colorbar;
caxis([0 1]);
set(hdl,'fontsize',14);
set(gca,'fontsize',14);
axis('square');

subplot(1,3,2);
imagesc(sphereAngs/pi,dtSamps,Pij);
xlabel('\bf{\phi/\Pi}','fontsize',14);
ylabel('\bf{t}','fontsize',14);
title('\bf{Projection data}','fontsize',14);
hdl = colorbar;
set(hdl,'fontsize',14);
set(gca,'fontsize',14);
axis('square');

subplot(1,3,3);
imagesc(Xmn(1,:),Ymn(:,1),real(Omn_r));
xlabel('\bf{x}','fontsize',14);
ylabel('\bf{y}','fontsize',14);
hold on;
text(Xmn(1,1),Ymn(6,1),sprintf('\\bf{N_{\\phi} = %d}',...
     nHalfSphereSamples),'fontsize',12,'color',[1 1 1]);
text(Xmn(1,1),Ymn(16,1),sprintf('\\bf{\\Delta_t = %.2f}',...
     sampDt),'fontsize',12,'color',[1 1 1]);
title('\bf{Reconstructed medium}','fontsize',14);
hold off;
hdl = colorbar;
caxis([0 1]);
set(hdl,'fontsize',14);
set(gca,'fontsize',14);
axis('square');

figHdl = figure;
imagesc(Xmn(1,:),Ymn(:,1),real(Omn_r));
xlabel('\bf{x}','fontsize',14);
ylabel('\bf{y}','fontsize',14);
hold on;
text(Xmn(1,1),Ymn(6,1),sprintf('\\bf{N_{\\phi} = %d}',...
     nHalfSphereSamples),'fontsize',12,'color',[1 1 1]);
text(Xmn(1,1),Ymn(16,1),sprintf('\\bf{\\Delta_t = %.2f}',...
     sampDt),'fontsize',12,'color',[1 1 1]);
title('\bf{Reconstructed medium}','fontsize',14);
hold off;
hdl = colorbar;
caxis([0 1]);
set(hdl,'fontsize',14);
set(gca,'fontsize',14);
axis('square');


print(sprintf('rec_dt_%.2f_Nphi_%d.eps',sampDt,nHalfSphereSamples),'-depsc');
      
close(figHdl);

