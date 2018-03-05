clc
clear

close all

more off

% Sampling distance along line. 
gMax = pi/8;
D_rad_rat = 1*sqrt(2)*1.25;
%D_rad_rat = 2;

nHalfSphereSamples = 150;

Omn = phantom(100);

Xrange = [-4,4];
Yrange = [-4,4];

nPts = size(Omn);

[Xmn,Ymn] = meshgrid(linspace(Xrange(1),Xrange(2),nPts(1)),linspace(Yrange(1),Yrange(2),nPts(2)));

%% Create set of projections

% Determine radius of object domain
R0 = sqrt((Xrange(2)-Xrange(1))^2+(Yrange(2)-Yrange(1))^2)/2;
D = D_rad_rat*R0;


nSphereSamples = nHalfSphereSamples*2;
dg = 2*pi/nSphereSamples;

[Rij,gSamps,D,b] = PWLequaAngRaysFrom2Dmesh_iterative(Xmn,Ymn,Omn,dg,gMax,D,nSphereSamples);

%Omn_r = unfilteredBackprojection_TD_equiAng(Rij,dg,gSamps,D,R0,b,Xmn,Ymn);
[Omn_r,Pij] = unfilteredBackprojection_TD_iradon(Rij,dg,gSamps,D,R0,b,Xmn,Ymn);


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
if(exist('Pij'))
  imagesc(b*180/pi,gSamps(:),Pij);
else
  imagesc(b*180/pi,gSamps(:),Rij);
end
xlabel('\beta^o');
ylabel('\gamma^o');
title('Projection data');
colorbar;
axis('square');

subplot(1,3,3);
imagesc(Xmn(1,:),Ymn(:,1),real(Omn_r));
xlabel('X');
ylabel('Y');
title('Reconstructed medium');
colorbar;
%caxis([0 1]);
axis('square');



