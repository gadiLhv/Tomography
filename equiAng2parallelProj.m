function [Pij,dt,tVals] = equiAng2parallelProj(Rij,gVals,D)

Ns = size(Rij,1);
% Number of samples for size
N = (Ns-1)/2;


ms = (0:2*N).';
ns = ms;

% Sampling indices vectors (mainly for filter creation)
n = (-N:N).';
k = n;

% Re-arrange samples
% Number of samples around circumference
M = size(Rij,2);
m = 0:M-1;



% phi sample number
p = mod(bsxfun(@plus,n,m),M);
countTs = sum(bsxfun(@eq,unique(p(:)).',p(:)),1);

% Validate that there are as many 't' values for each phi sample
if(min(countTs) ~= max(countTs)) error('Unable to find a valid number of samples'); end

countTs = countTs(1);
numCols = numel(Rij)/countTs;

resortIdxs = bsxfun(@plus,p*countTs,(1:countTs).');
resortIdxs = resortIdxs(:);

Pij = reshape(Rij(resortIdxs),...
              [countTs numel(Rij)/countTs]);

% Now, since the sampling is not uniform along the t-axis,
% interpolate the necessary values in uniform sampling spaces.
tVals_r = D*sin(gVals);
dt = (max(tVals_r)-min(tVals_r))/(numel(gVals)-1);

% Create t value array to interpolate
tVals = linspace(min(tVals_r),max(tVals_r),numel(gVals));

% Array specific for interpolation
Pij_I = zeros([numel(tVals) size(Pij,2)]);

% Interpolate each column seperately
for pIdx = 1:size(Pij,2)
  Pij_I(:,pIdx) = interp1(tVals_r,...
                          Pij(:,pIdx),...
                          tVals(:),...
                          'cubic',0);
end

Pij = Pij_I;

end