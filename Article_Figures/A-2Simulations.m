function [ Err] =figure4Simulations(noiseLevel,sampleLocation,PDmodel,FitModels,AddFitPolyorder,iter,lambdaRidge,lambdaT1,lambdaRidgeBL) 
% the function simulate the data and estimate the PD with different
% regularization method. and calculate the goodness of estimation.
%% Supplementary Figure
%
% Compare PD estimation using different types of regression constraints
%
% We compare ridge, T1 and coil correlation and no regularization at all.
%
% The conclusion is that T1-regularization has the best overall properties.
% The absolute error is slightly smaller with ridge regression, but the
% spatial distribution of the error is systematic over the volume with
% ridge, and the spatial error is white noise with T1 regularization.
%
% No regularization is poor.  Coil correlation is OK, but not quite as good
% in some instances as T1 regularization.
%
% FitModels =1 a vector 6x1. one to use a model zero for not
% FitModels= [Cor        Ridge    Cor+Ridge     T1- reg       No-Reg       T1-Re+Tissue-mask]
%FitModels= [1 1 1 1 1 0] (defult)
%
% Requires vistasoft, mrQ and knkutils
%
% AM/BW Vistaosft Team, 2013

%%  Make sure mrQ is on the path
%
% addpath(genpath(fullfile(mrqRootPath)));
% addpath(genpath('/Users/wandell/Github/knkutils')); gitRemovePath;

%% Generate example parameters for the coils from the phantom data

if notDefined('FitModels')
    FitModels=ones(6,1);
    FitModels(6)=0;
end

if notDefined('AddFitPolyorder')
AddFitPolyorder=1;
end

if notDefined('iter')
iter=1000;
end

nSamples = 3;      % The box is -nSamples:nSamples
nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
pOrder   = 2;      % Second order is good for up to 5 samples
noiseFloor = 500;  % This is the smallest level we consider
%sampleLocation = 3;% Which box
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data
BasisFlag    = 'qr';    % Which matrix decomposition for fitting.

% This produces the key parameters for the polynomial approximations.  The returned variables includes
% the polynomial basis, pBasis, the M0 data, M0S_v, additional parameters,
% such as the box size.
phantomP = pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
    noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);

[pBasis]  = polyCreateMatrix(nSamples,pOrder+AddFitPolyorder,nDims,BasisFlag);

%% Simulate coil gain using the poylnomial fits to the phantom data

 nUseCoils = 4;                         % How many coils to use
MaxcoilNum=16;                     %last coil to consider
% These are typical coil functions
% We can sort coils by minimalcorrelation between the coils to find the best set.
% We use this algorithm to select the coils
% Find the minimum correlation (min abs corr give us the set with corr that
% are closer to zero).  Choose those coils.

coils=mrQ_select_coilsMinCorrelation(nUseCoils,MaxcoilNum,phantomP.M0_v);

% Get the poylnomial coeficents for those coils
GainPolyPar = phantomP.params(:,coils);

% Create the coil gains over voxels by multiplying the polynomial
% coeficents and the polynomial basis.
G = phantomP.pBasis*GainPolyPar;

%% Simulate PD

% There are several PD spatial types.
% Type help mrQ_simulate_PD
[PD, R1, Tissuemask] = mrQ_simulate_PD(PDmodel,phantomP.nVoxels);
Err.SimPD=PD(:);
Err.SimR1=R1;
Err.Tissuemask=Tissuemask(:);

 %showMontage(PD)

%% Simulate MRI SPGR signal with noise

%noiseLevel = 2;   % ?? Units???

% Simulate the M0 and T1 fits of multi SPGR images.
[MR_Sim] = simSPGRs(G,PD(:),[],R1(:),[],[],noiseLevel,printImages);
Err.R1Fit=MR_Sim.R1Fit;
% MR_Sim is a structure with multiple fields that include the simulation
% inputs MR sigunal inputs and the calculations from fitting the signal
% equation.

%% Initiate fit params
kFold   = 2; % X-validate on half the data

% The nonlinsqr fit with ridge fits
[PDinit, g0]=Get_PDinit(0,[],1,MR_Sim.M0SN,pBasis);

options = optimset('Display','off',...
    'MaxFunEvals',Inf,...
    'MaxIter',iter,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','levenberg-marquardt');
boxSize = repmat(phantomP.rSize,1,nDims);

%[Poly1,str] = constructpolynomialmatrix3d(boxSize,find(ones(boxSize)),1);



%% Coil correlation regularization
if FitModels(1)==1;
coefdat = tril(corrcoef(MR_Sim.M0SN),-1);

RegWeight  = 1000;
TissueMask = logical(MR_Sim.M0SN(:,1));
XvalidationMask = logical(MR_Sim.M0SN);

[g, resnorm,dd1,exitflag] = ...
    lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,MR_Sim.M0SN,pBasis,nUseCoils,RegWeight,TissueMask,coefdat,XvalidationMask),g0,[],[],options);

[PD_cor, Gn] = pdEstimate(MR_Sim.M0SN,pBasis, g);

scale     = mean(PD(:)./PD_cor(:));
PD_cor    = PD_cor(:)*scale;
Err.Err_cor(1)    = (calccod(PD_cor(:),PD(:)));  %R^2
Err.Err_cor(2)    = 100*mean(abs(PD(:)-PD_cor(:))./(PD(:))); % mean abs err
Err.Err_cor_G=g;

%figure;plot(PD(:),PD_cor(:),'*');identityLine(gca);axis image; axis square
PD_cor    = reshape(PD_cor,boxSize);

%showMontage(PD_cor)

ErrMap=PD-PD_cor;
[~,~,Err.Err_cor_LinErrResid] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);

end




%% Ridge regularization X- validation
if FitModels(2)==1;
    
%  To make
end
%% Ridge & Bi-Linear X- validation
if FitModels(3)==1;
   
maxLoops = 1000;
sCriterion = 1e-4;  % Stopping criterion    

  if notDefined('lambdaRidgeBL')
        lambda = [10 5 2 1 0.1 0.05 0.01 0.005 0.001 0];
    else
        lambda =lambdaRidge;
        
  end
    
  if  length(lambda)>1
% Set the fiiting loop parmeters
  [X_valdationErr,  ]=pdX_valdationLoop_RidgeReg_ver1(lambda,kFold,MR_Sim.M0SN,pBasis,PDinit,maxLoops,sCriterion);
  
  %mrvNewGraphWin;loglog(lambda+1e-1,X_valdationErr(2,:),'*-'); xlabel('lambda');ylabel('X-V error');

BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)));
RidgeBLFit = pdBiLinearRidgeFit(MR_Sim.M0SN,pBasis,lambda(BestReg),maxLoops,sCriterion,PDinit);
Err.Err_ridgeBL_XV_Lamda=lambda(BestReg);

  else
    RidgeBLFit = pdBiLinearRidgeFit(MR_Sim.M0SN,pBasis,lambda,maxLoops,sCriterion,PDinit);

end

  % Use the best lambda and fit the full data set

scale      = mean(PD(:)./RidgeBLFit.PD(:));
PD_ridgeBL_XV   = RidgeBLFit.PD(:)*scale;
%PD_ridgeBL_XV     = reshape(PD_ridgeBL_XV,boxSize);
%showMontage(PD_ridgeBL_XV)
%figure;plot(PD(:),PD_ridgeBL_XV(:),'*');identityLine(gca);axis image; axis square

Err.Err_ridgeBL_XV (1)    = (calccod(PD_ridgeBL_XV(:),PD(:))); %R^2
Err.Err_ridgeBL_XV (2)    =  100*mean(abs(PD(:)-PD_ridgeBL_XV (:))./(PD(:))); % mean abs err
Err.Err_ridgeBL_XV_G=RidgeBLFit.g;

  PD_ridgeBL_XV   = reshape(PD_ridgeBL_XV,boxSize);

ErrMap=PD-PD_ridgeBL_XV;
[~,~,Err.Err_ridgeBL_LinErrResid] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);


end
%% R1 regularization
if FitModels(4)==1;
    if notDefined('lambdaT1')
        % Possible weights to test.  We will choose the one that cross-validates best.
        lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0];    else
        lambda =lambdaT1;
        
    end
    
    

Rmatrix(1:phantomP.nVoxels,1) = 1;
Rmatrix(:,2) = double(MR_Sim.R1Fit);

if  length(lambda)>1

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,pBasis,Rmatrix,g0,[],[]);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)));

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg),[],options);
Err.Err_T1reg_Lamda=lambda(BestReg);

else
    
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda,MR_Sim.M0SN,pBasis, ...
    Rmatrix, g0,[],options);

end


scale     = mean(PD(:)./NL_T1reg.PD(:));
PD_T1reg  = NL_T1reg.PD(:)*scale;
%PD_T1reg  = reshape(PD_T1reg,boxSize);
%showMontage(PD_T1reg)

Err.Err_T1reg(1)    = (calccod(PD_T1reg(:),PD(:))); %R^2
Err.Err_T1reg(2)    =  100*mean(abs(PD(:)-PD_T1reg(:))./(PD(:))); % mean abs err
%figure;plot(PD(:),PD_T1reg(:),'*');identityLine(gca);axis image; axis square

Err.Err_T1reg_G=NL_T1reg.g;


PD_T1reg  = reshape(PD_T1reg,boxSize);

ErrMap=PD-PD_T1reg;
[~,~,Err.Err_T1reg_ErrResid] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);

Err.PD_T1reg=PD_T1reg;

end
%% no regularization
if FitModels(5)==1;
NL_NoReg  = pdBiLinearFit_lsqSeach(MR_Sim.M0SN,pBasis,PDinit,options);
scale     = mean(PD(:)./NL_NoReg.PD(:));
PD_Noreg  = NL_NoReg.PD(:)*scale;
% showMontage(PD_Noreg)
Err.Err_Noreg(1)    = (calccod(PD_Noreg(:),PD(:))); %R^2
Err.Err_Noreg(2)    =  100*mean(abs(PD(:)-PD_Noreg(:))./(PD(:))); % mean abs err

PD_Noreg  = reshape(PD_Noreg,boxSize);
ErrMap=PD-PD_Noreg;
[~,~,Err.Err_Noreg_ErrResid] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);

Err.PD_Noreg=PD_Noreg;

Err.Err_Noreg_G=NL_NoReg.g;

%figure;plot(PD(:),PD_Noreg(:),'*');identityLine(gca);axis image; axis square
end
%% R1 regularization + tissue mask
if FitModels(6)==1;
% Possible weights to test.  We will choose the one that cross-validates best.
lambda = [1e4 5e3 1e3 5e2 1e2 5e1 1e1 5e0 1e0 5e-1 1e-1 0];
Rmatrix(1:phantomP.nVoxels,1) = 1;
Rmatrix(:,2) = double(MR_Sim.R1Fit);

% Loop over regularization weights and calculate the X-validation error
[X_valdationErr,   gEstT ] = ...
    pdX_valdationLoop_2(lambda,kFold,MR_Sim.M0SN,pBasis,Rmatrix,g0,Tissuemask,options);

% Find the lambda that best X-validates (minimal RMSE error)
BestReg = find(X_valdationErr(2,:) == min(X_valdationErr(2,:)));

% Use the best lambda and fit the full data set
[NL_T1reg.PD,~,NL_T1reg.G,NL_T1reg.g, NL_T1reg.resnorm,NL_T1reg.exitflag ] = ...
    pdCoilSearch_T1reg(lambda(BestReg),MR_Sim.M0SN,pBasis, ...
    Rmatrix, gEstT(:,:,1,BestReg),Tissuemask,options);

scale     = mean(PD(:)./NL_T1reg.PD(:));
PD_T1reg_TM  = NL_T1reg.PD(:)*scale;
PD_T1reg_TM  = reshape(PD_T1reg_TM,boxSize);
%showMontage(PD_T1reg)

Err.Err_T1reg_TM  (1)    = (calccod(PD_T1reg_TM  (:),PD(:))); %R^2
Err.Err_T1reg_TM  (2)    =  100*mean(abs(PD(:)-PD_T1reg_TM  (:))./(PD(:))); % mean abs err

ErrMap=PD-PD_T1reg_TM;
[~,~,Err.Err_T1reg_TM_ErrResid] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);

Err.PD_T1reg_TM=PD_T1reg_TM;

Err.Err_T1reg_TM_G=NL_T1reg.g;


%figure;plot(PD(:),PD_T1reg(:),'*');identityLine(gca);axis image; axis square
end

%% Notes

%% END


return
% close all
% PD_T1reg  = reshape(PD_T1reg,boxSize);
% [ErrFFT_R1,autoCor_R1]=CalBias_Noise(PD,PD_T1reg);
% showMontage(ErrFFT_R1);title('FFT R1 reg');showMontage(autoCor_R1);title('Auto Cor R1 reg')
% 
% PD_ridgeBL_XV=reshape(PD_ridgeBL_XV,boxSize);
% [ErrFFT_Rd,autoCor_Rd]=CalBias_Noise(PD,PD_ridgeBL_XV);
% showMontage(ErrFFT_Rd);title('FFT R1 ridge');showMontage(autoCor_Rd);title('Auto Cor  ridge')
