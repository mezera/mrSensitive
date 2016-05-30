function  [X_valdationErr,   BLFit_RidgeReg, FitT, useX, kFold,lambda1, maxLoops,sCriterion]=pdX_valdationLoop_RidgeReg_ver1(lambda1,kFold,M0,pBasis,PDinit,maxLoops,sCriterion)
% Fit M0 for coil gain and PD with PD with different weight (lambda1) for PD
% ridge by T1 (linear relations).
%
%  [X_valdationErr   gEstT, resnorm, FitT]= ...
%    pdX_valdationLoop_RidgeReg_ver1( lambda1,kFold,M0,pBasis,GainPolyPar,maxLoops,sCriterion)
%
% Calculates the Cross Validation eror for each regularization wight
%
% AM  & BW VISTASOFT Team, 2013


%% intilaized parameters
nVoxels=size(M0,1);
Ncoils=size(M0,2);
nPolyCoef=size(pBasis,2);

if notDefined('PDinit')
[PDinit]=Get_PDinit(0,[],1,M0,pBasis);
end
if notDefined('kFold')
kFold=2;
end
if notDefined('lambda1')
lambda1=[10 5 2 1 0.1 0.05 0.01 0.005 0.001 0];
end

% Stopping criterion  

if notDefined('maxLoops')
maxLoops = 1000;
end
if notDefined('sCriterion')
sCriterion = 1e-4;  
end

%% intilaized X-Validation and fits

% separate the data to Kfold fit and X-Validation part
[ useX, kFold] = getKfooldCVvoxel_full(nVoxels,Ncoils,kFold);

% argumante to save the X-Validation results
X_valdationErrKfold = zeros(2,kFold);
X_valdationErr  = zeros(2,length(lambda1));


%% loop over lambda1 regularization wight
% Sweep out the lambda values
for ii=1:length(lambda1),
    
    % Loop over the kFold Cross Validation
    for jj=1:kFold
        %select the position to estimate the function
       
        % Searching on the gain parameters, G.
        
                FitMask=zeros(size(M0));FitMask(find(useX~=jj))=1;FitMask=logical(FitMask);
  
       BLFit_RidgeReg(jj,ii) = pdBiLinearRidgeFit_XV(M0,pBasis,lambda1(ii),maxLoops,sCriterion,PDinit,FitMask);
           
       
        
         %%  calculate X-Validation error:
        
        Xmask=zeros(size(M0));Xmask(find(useX==jj))=1;Xmask=logical(Xmask); % Xmask are the hold position form X-Validation.
        
                %Check if the coil coefficent can explain the hold data
        [FitT(jj,ii).err_X, FitT(jj,ii).err_F] = X_validation_errHoldVoxel_full( BLFit_RidgeReg(jj,ii).g ,M0,pBasis,nVoxels,Ncoils,FitMask,Xmask);
        
        % Two posible error function. we don't find a big different between
        % them
        X_valdationErrKfold(1,jj)= sum(abs( FitT(jj,ii).err_X (:))); %sum of absulot erro
        X_valdationErrKfold(2,jj)= sum(  FitT(jj,ii).err_X (:).^2); %RMSE
    end
    
    %sum over the Kfold X-Validation error for this lambda1
    X_valdationErr(1,ii)=sum(X_valdationErrKfold(1,:));
    X_valdationErr(2,ii)=sum(X_valdationErrKfold(2,:));
    
    
end

