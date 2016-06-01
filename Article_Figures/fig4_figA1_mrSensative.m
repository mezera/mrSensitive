% Figure 4 & Apendix figure 1-A 
%
% Illustrate how PD can be measure from M0 while estimate the coil gain using
%  bi-linear solutions
%
% The problem we pace is the effect of noise on this solutions
%
% % AM/BW  Mezer Lab Wandell Lab & Vistaosft Team, 2013

%% 
clear
disp('Figure 4 and Apendix figure 1-A ')

%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));

%% Generate example parameters for the coils from the phantom data

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
sampleLocation = 4;% Which box
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations.  The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
boxSize = repmat(phantomP.rSize,1,nDims);

%% Simulate the proton density (PD)
[PD, R1] = mrQ_simulate_PD('6',phantomP.nVoxels);


%% Simulate coil gain 
% We use the poylnomial fits to the phantom data a typical coil function

% Select a set of coils
coils = [1 3 5 8];
% Get those coil poylnomyal coeficents
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomials
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;

%% Simultae MRI SPGR signal and fit with and without noise
noiseLevel = 2;   % ?? Units???

% Simultate the M0 and T1 fits of multi SPGR images.
% We should reset the random noise generator here to the same value before
% making this call.
rng('default')
[MR_Sim]= simSPGRs(G,PD(:),[],[],[],[],noiseLevel,false);

%% Separate the M0 into (PD,Gain) using bilinear alternating least squares (ALS) or non linear search

% No noise
lambda=0; % this is the general call for ridge lsq search. we use it here with no regularization.

% it is the non linear search. it is much faster than the Alternative least square  (see below)
BLFit_Fit_NoNoise = pdBiLinearFit_lsqRidgeSeach(MR_Sim.M0S,phantomP.pBasis,[],[],lambda);

scale      = mean(PD(:)./BLFit_Fit_NoNoise.PD(:));
PD_NoNoise   = BLFit_Fit_NoNoise.PD(:)*scale;
%% similar result can be reached with ALS but it will be much much longer
% we comment this as it will take hours. 
%
%     maxloop=1.5e6; % many iterations are needed
%     RidgeBLFit = pdBiLinearRidgeFit(MR_Sim.M0S,phantomP.pBasis,lambda,maxloop,1e-8); %ALS (ridge) with Lamda of 0 no regularization.
%    scale      = mean(PD(:)./RidgeBLFit.PD(:)); 
%    PD_NoNoiseBL   = RidgeBLFit.PD(:)*scale;



%% Scatter plot of  estimated PD without noise vs. simulated PD

MM = minmax([ PD(:) PD_NoNoise(:)]);
mrvNewGraphWin('Figure 4','wide');
subplot(1,3,1)
hold on
plot(PD_NoNoise(:),PD(:),'o' ,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')

%if using the ALS
%plot(PD_NoNoiseBL(:),PD(:),'o' ,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')

xlabel('Estimated PD','FontSize',16); ylabel('True PD','FontSize',16);
xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)])
axis image; axis square
%legend('PD estimate without noise','Location','NorthWest')
title ('figure 4a: Noise free','FontSize',16)
xlim([0 2]);ylim([0 2])
identityLine(gca);
set(gca,'FontSize',16)

%% With noise
lambda=0; % this is the general cool for ridge regretion. we use it here with no regularization.

BLFit_Fit_Noise = pdBiLinearFit_lsqRidgeSeach(MR_Sim.M0SN,phantomP.pBasis,[],[],lambda);

scale      = mean(PD(:)./BLFit_Fit_Noise.PD(:));
PD_Noise   = BLFit_Fit_Noise.PD(:)*scale;

%%  Scatter plot of PD estimate with  noise vs Simulated PD


MM = minmax([ PD(:) PD_Noise(:)]);
%mrvNewGraphWin;
subplot(1,3,2)

hold on
plot(PD_Noise(:),PD(:),'o' ,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')

xlabel('Estimated PD','FontSize',16); ylabel('True PD','FontSize',16);
xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)])
axis image; axis square
%legend('PD estimate with noise','Location','NorthWest')
title ('figure 4b: Noise overfit','FontSize',16)
xlim([0 2]);ylim([0 2])
identityLine(gca);
set(gca,'FontSize',16)

%% T1 regularization
% make the T1 regularization basis 
Rmatrix(1:phantomP.nVoxels,1) = 1;
Rmatrix(:,2) = double(MR_Sim.R1Fit);
%X- validate the best regularization value
kFold=2;
lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0];
[X_valdationErr,   gEstT ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,phantomP.pBasis,Rmatrix,[],[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)));

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,phantomP.pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg),[],[]);

scale      = mean(PD(:)./NL_T1reg.PD(:));
PD_T1   = NL_T1reg.PD(:)*scale;


%%  Scatter plot of PD estimate with  noise and T1 regularization vs Simulated PD

MM = minmax([ PD(:) PD_T1(:)]);
%mrvNewGraphWin;
subplot(1,3,3)

hold on
plot(PD_T1(:),PD(:),'o' ,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')

xlabel('Estimated PD','FontSize',16); ylabel('True PD','FontSize',16);
xlim([MM(1) MM(2)]); ylim([MM(1) MM(2)])
axis image; axis square
%legend('PD estimate with noise T1 reg','Location','NorthWest')
xlim([0 2]);ylim([0 2])
identityLine(gca);
set(gca,'FontSize',16)
title ('figure 4c :T1 regularization','FontSize',16)


%% Validation error

disp( 'Figure A-1')

% log log
% mrvNewGraphWin;
% %
% lambda(end)=lambda(end)+1e-10;
% loglog(lambda,X_valdationErr(2,:),'-ko' ,'MarkerSize',5,'MarkerFaceColor','k'); 
% xlabel('regularization weight ( \lambda)');ylabel('Cross validation (SSE) '); 
% set(gca,'ylim',[10^5.3 10^6.5])
% %set(gca,'YTick',[10^5.4  ;10^6 ;10^6.4])
% %set(gca,'YTickLabel',{'10 ^5.4'  ;'10 ^6' ;'10 ^6.4'})
% %set(gca,'XTick',[0 1000 3000 5000 7000 9000])
% %set(gca,'XTickLabel',{'0' '1000' '3000' '5000' '7000' '9000'})
% 
% %set(gca,'xlim',[-500 10400])
% 
% set(gca,'xlim',[10^-1.2 10^4.2])
% grid on
% % Set yticks now ...
% % set(gca,'ytick',[10^5.5 10^6 10^6.5])



%semilog
mrvNewGraphWin('Apendix Figure 1-A');
%
semilogy(lambda,X_valdationErr(2,:),'-ko' ,'MarkerSize',8,'MarkerFaceColor','w','linewidth',2); 
xlabel('regularization weight ( \lambda)','FontSize',20);   ylabel('Cross validation (SSE) ','FontSize',20); 
set(gca,'XTick',[0 1000 3000 5000 7000 9000])
set(gca,'XTickLabel',{'0' '1000' '3000' '5000' '7000' '9000'})
set(gca,'FontSize',20)

set(gca,'xlim',[-500 10400])

grid on
set(gca,'ylim',[10^5.5 10^6.5])



%set(gca,'YTick',[  10^5.5  10^5.6 10^5.7 10^5.8 10^5.9 10^6  10^6.1  10^6.2  10^6.3   10^6.4])
%set(gca,'YTickLabel',{'' '' '10 ^5.7'  '' '' '' '' '10 ^6' '' ''  '' ''})
