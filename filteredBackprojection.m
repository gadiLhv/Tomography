function Omn_r = filteredBackprojection(Pij,sampDt,filterBW,tVals,phiVals,Xmn,Ymn)

% Create slices by performing a Fourier transform

% The formula for FFT is F = sum_j(f(j)*exp(-1i*2*pi/N)^(j*k))
% Corresponding IFFT is f = 1/N*sum_k(F(k)*exp(1i*2*pi/N)^(j*k))
Ns = (size(Pij,1)-1)/2;
Ts = 1/sampDt;

ms = (0:2*Ns).';
ns = ms;


P_shifted = bsxfun(@times,Pij,exp(1i*2*pi*Ns*ns(:)/(2*Ns+1)));
% Slices:
S = bsxfun(@times,...
           exp(-1i*2*pi*Ns*(Ns-ms(:))/(2*Ns+1)),...
           fft(P_shifted,2*Ns+1,1));

     
S_shifted = bsxfun( @times,...
                    S,...
                    exp(-1i*2*pi*Ns*ms(:)/(2*Ns+1)));

% Unfiltered slice quantity
S_shifted_UF =  bsxfun( @times,...
                        S_shifted,...
                        (Ts/(2*Ns+1))*abs(ms(:)-Ns));



                        
% Number of frequency samples composing the wind
Nw = ceil(Ns*filterBW);
% Index of middle frequency
Nmid = (numel(tVals)-1)/2+1;

%%% Construct the hanning window
%%Hwin = zeros([numel(tVals) 1]);
%%Hwin(Nmid-Nw:Nmid+Nw) = 1;
%%Hwin = Hwin(:);
%%
%%
%%                
%%% Filter Fourier slices
%%S_shifted_Filtered = bsxfun(@times,S_shifted_UF,Hwin);
%%
%%
%%figure;
%%imagesc(abs(S_shifted_Filtered));

% The reconstructed quantity
Q = bsxfun(@times,...
           exp(1i*2*pi*Ns*(Ns-ns(:))/(2*Ns+1)),...
           ifft(S_shifted,2*Ns+1,1));


           
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