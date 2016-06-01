%% Figure 1
%
% Illustrate how PD and coil sensitivity are joined to make the M0  images %
% Code associated with Mezer, et. al. 2016, HBM
% 
%
% % AM/BW  Mezer Lab Wandell Lab & Vistaosft Team, 2013


%% 
disp('Figure 1')

%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));

%% Generate example parameters for the coils from the phantom data

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 2;% Which box
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations.  The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
boxSize = repmat(phantomP.rSize,1,nDims);

%% simulate PD
[PD, R1] = mrQ_simulate_PD('6',phantomP.nVoxels);


%% Simulate coil gain 
% We use the poylnomial fits to the phantom data a typical coil function

% Select a set of coils
coils = [1 3 5 ];
% Get those coil poylnomyal coeficents
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomials
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;



%% Simultae MRI SPGR signal and fit without noise
noiseLevel = 0;   % ?? Units???

% Simultate the M0 and T1 fits of multi SPGR images.
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,false);

%% Show an example slice of PD the map

mrvNewGraphWin;

slice = 4;
imagesc(PD(:,:,slice));
colormap(gray); axis image; axis off
title('Simulated PD');
%% The M0 images
mrvNewGraphWin([],'tall');

mn = min(MR_Sim.M0S(:)); mx = max(MR_Sim.M0S(:));
for ii=1: length(coils)
    subplot(4,1,ii)
    M0=MR_Sim.M0SN(:,ii);
    M0=reshape(M0,boxSize);
    imagesc(M0(:,:,slice));
    caxis([mn mx]);
    colormap(gray); axis image; axis off;
    title(sprintf('M0 for coil %d\n',ii));
    
    % mrUtilResizeFigure(gcf, 900, 900);
    % mrUtilPrintFigure(['M0_example_slice' num2str(ii) '.eps']);
end



%% The G images
mrvNewGraphWin([],'tall');

mn = min(G(:)); mx = max(G(:));
for ii=1: length(coils)
    subplot(4,1,ii)
    Gi=G(:,ii);
    Gi=reshape(Gi,boxSize);
    imagesc(Gi(:,:,slice));
    caxis([mn mx]);
    colormap(gray); axis image; axis off;
    title(sprintf('Sensativity for coil %d\n',ii));
    
    % mrUtilResizeFigure(gcf, 900, 900);
    % mrUtilPrintFigure(['M0_example_slice' num2str(ii) '.eps']);
end
%% End
