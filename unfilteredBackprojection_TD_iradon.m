function [Omn_r,Pij] = unfilteredBackprojection_TD_iradon(Rij,dg,gVals,D,R0,bVals,Xmn,Ymn)

[Pij,dt,tVals] = equiAng2parallelProj(Rij,gVals,D);

Ns = size(Pij,1);
% Number of samples for size
N = (Ns-1)/2;


ms = (0:2*N).';
ns = ms;

% Sampling indices vectors (mainly for filter creation)
n = (-N:N).';
k = n;

% Re-arrange samples
% Number of samples around circumference
M = size(Pij,2);
m = 0:M-1;

% From this point on this goes like a regular parallel-line projection
tau = dt;
kappa = 1/tau;

% Create filter
Hn = (mod((-N:N).',2) ~= 0).*(-1./((pi*tau*(-N:N).').^2));
Hn(N+1) = (1/(4*tau^2));

% Re-organize the filter in a cyclic fashion (instead of multiplying by exponents)
Hn = [Hn(N+1:end) ; Hn(1:N)];
% Do the same with the projections
Pij = [Pij(N+1:end,:) ; Pij(1:N,:)];

% Pad with zeros if necessary (this is future compatibility with additional filters)
maxSize = max([size(Hn,1) ; size(Pij,1)]);
Npad_Hn = max([maxSize-size(Hn,1)-1 ; 0]);
Npad_Pij = max([maxSize-size(Pij,1)-1 ; 0]);

Hn = [Hn ; zeros([Npad_Hn 1])];
Pij = [Pij ; zeros([Npad_Pij size(Pij,2)])];

% Fourier transform the padded filter
Hk = fft(Hn,size(Hn,1),1);

% Fourier transform the padded projections (into padded slices)
Sij_k = fft(Pij,size(Pij,1),1);

% Multiply and shift back in time:
Qij_k = bsxfun(@times,Sij_k,Hk);

% Create hamming window and shift it
K = hamming(size(Qij_k,1));
Nshift = floor(size(Qij_k,1)/2)+1;
K = [K(end-Nshift+1:end) ; K(1:Nshift-1)];
                          
Q = tau*ifft(Qij_k,size(Hn,1),1);

Q = Q(1:Ns,:);

% Re organize the reconstruction quantity
Q = [Q(N+2:end,:) ; Q(1:N+1,:)];
Pij = [Pij(N+2:end,:) ; Pij(1:N+1,:)];

phiVals = bVals;

% Now it is time to interpolate the (x,y) values of the object domain
tVals = tVals(:);       % Used to interpolate the interpolated t values 
phiVals = phiVals(:).';     % Used for summation


% The required interpolated t values
tI = bsxfun(@times,Xmn(:),cos(phiVals+pi)) + bsxfun(@times,Ymn(:),sin(phiVals+pi));



% Interpolate Q values, accordingly
QI = zeros(size(tI));
for angIdx = 1:size(tI,2)
  QI(:,angIdx) = interp1(tVals,...
                         Q(:,angIdx),...
                         tI(:,angIdx),...
                         'Cubic',0);
end


Omn_r = (pi/numel(phiVals))*sum(QI,2);

Omn_r = reshape(Omn_r,size(Xmn));