function Omn_r = unfilteredBackprojection(Pij,Dt,tVals,phiVals,Xmn,Ymn)
% Omn_r = unfilteredBackprojection(Pij,Dt,tVals,phiVals,Xmn,Ymn)
%
% Performs a filtered backprojection discrete approximation.
% 
% Inputs: 
% Pij - [M,N] complex value matrix denoting the projections. Along the 
%       columns varies the sample along the 't' direction (tangential)
%       while along the columns the projectioin direction varies.
% Dt -  Sampling rate. Scalar.
% tVals - [M,1] real valued vector denoting the samples along the 't'
%       direction.
% phiVals - [N,1] real valued vector denoting the directions in which
%       the projections are made, between [0,pi).
% Xmn - [Ms,Ns] real valued matrix of 'x' coordinates in which the 
%       reconstruction needs to be calculated.
% Ymn - Same as 'Xmn' but with 'y' coordinates.
%
% Outputs:
% Omn_r - [Ms,Ns] matrix of the contrast profile reconstruction. 

% Calculate "half" the number of samples
M = (size(Pij,1)-1)/2;

% Calculate the maximal k (sampling frequency)
k_M = 1/(2*Dt);

% Explicitly write sampling rates in both frequency and 
% spatial domain
Dk = (2*k_M)/(2*M + 1);

% Fourier transform matrix (denoted as E-)
E = exp(-1i*2*pi*(-M:M).'*(-M:M)/(2*M+1));

% Slices (sinograms):
S_UF = Dt*E*Pij;

% Filter the slices
S_F = bsxfun(@times,S_UF,abs(Dk*(-M:M).'));

% Reconstruction quantity is obtained with the IFT matrix
Q = Dk*E'*S_F;



% Now it is time to interpolate the (x,y) values of the object domain
tVals = tVals(:);       % Used to interpolate the interpolated t values 
phiVals = phiVals(:).';     % Used for summation

% The required interpolated t values
tI = bsxfun(@times,Xmn(:),cos(phiVals+pi)) + bsxfun(@times,Ymn(:),sin(phiVals+pi));

% Interpolate Q values, accordingly
QI = zeros(size(tI));
for angIdx = 1:size(tI,2)
  QI(:,angIdx) = interp1(tVals,Q(:,angIdx),tI(:,angIdx),'linear',0);
end


Omn_r = (pi/numel(phiVals))*sum(QI,2);

Omn_r = reshape(Omn_r,size(Xmn));