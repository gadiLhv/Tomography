function Omn_r = unfilteredBackprojection_TD_equiAng(Rij,dg,gVals,D,R0,bVals,Xmn,Ymn)

Ns = size(Rij,1);
% Number of samples for size
N = (Ns-1)/2;

ms = (0:2*N).';
ns = ms;

% Sampling indices vectors (mainly for filter creation)
n = (-N:N).';
k = n;

% Create filter
gn = (mod((-N:N).',2) ~= 0).*((1./(pi*sin(gVals(:)))).^2);
gn(N+1) = (1/(8*dg^2));

% "Weight" the filter with the jacobian
Rij = D*bsxfun(@times,Rij,cos(gVals(:)));

% Re-organize the filter in a cyclic fashion (instead of multiplying by exponents)
gn = [gn(N+1:end) ; gn(1:N)];

% Do the same with the projections
Rij = [Rij(N+1:end,:) ; Rij(1:N,:)];

% Pad with zeros if necessary (this is future compatibility with additional filters)
maxSize = max([size(gn,1) ; size(Rij,1)]);
Npad_gn = max([maxSize-size(gn,1)-1 ; 0]);
Npad_Rij = max([maxSize-size(Rij,1)-1 ; 0]);

gn = [gn ; zeros([Npad_gn 1])];
Rij = [Rij ; zeros([Npad_Rij size(Rij,2)])];

% Fourier transform the padded filter
gk = fft(gn,size(gn,1),1);

% Fourier transform the padded projections (into padded slices)
Sij_k = fft(Rij,size(Rij,1),1);



% Multiply (convolve in the frequency domain)
Qij_k = bsxfun(@times,Sij_k,dg*gk);

figure;
subplot(3,2,1),imagesc(20*log10(abs(Sij_k))),colorbar,axis('square');
subplot(3,2,2),imagesc(20*log10(abs(Qij_k))),colorbar,axis('square');
subplot(3,1,2),plot(abs(gk));
subplot(3,1,3),plot(gn);

% Inverse Fourier transform to return to TD representation
Q = ifft(Qij_k,size(gn,1),1);

% Take only relevant (non-padded) part of the samples
Q = Q(1:Ns,:);

% Re organize the reconstruction quantity
Q = [Q(N+2:end,:) ; Q(1:N+1,:)];
           
% Now it is time to interpolate the (x,y) values of the object domain
gVals = gVals(:);       % Used to interpolate the interpolated t values 
bVals = bVals(:).';     % Used for summation

% In order to calculate the angles relative to the beam initial, they first have
% to be rotated to there

% The required interpolated gamma values
tI = bsxfun(@times,Xmn(:), cos(bVals+pi)) + bsxfun(@times,Ymn(:),sin(bVals+pi));
sI = bsxfun(@times,Xmn(:),-sin(bVals+pi)) + bsxfun(@times,Ymn(:),cos(bVals+pi)) + D;
gI = pi/2 - atan2(sI,tI);
LI2 = tI.^2 + sI.^2;

% Interpolate Q values, accordingly
QI = zeros(size(gI));
for angIdx = 1:size(gI,2)
  QI(:,angIdx) = interp1(gVals,...
                         Q(:,angIdx),...
                         gI(:,angIdx),...
                         'Cubic',0);
end

Omn_r = (2*pi/numel(bVals))*sum(QI./LI2,2);

Omn_r = reshape(Omn_r,size(Xmn));