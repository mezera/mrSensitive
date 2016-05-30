function err = errFitRidgeCorrNestBiLinear_lsq_v1(g,M0,pBasis,nCoils,W,RegWeight,coefdat,M0mask)
%Error function for bilinear coil gain estimate for lsq search with ridge
%error term and coil corralation term.
%
%  err = errFitRidgeCorrNestBiLinear_lsq(g,M0,pBasis,nVoxels,nCoils,W,RegWeight,coefdat)
%
% g       - Coil gain polynomial coefficients (nCoeffs x nCoils)
% M0      - Measurement (nVoxels x nCoils)
% pBasis  - Polynomial basis in columns (nVoxels x nCoeffs)
% nVoxels - Number of voxels
% nCoils  - Number of measurement coils
% W       - The ridge regression coefficients
% RegWeight The wight on the corralation panelty
% coefdat the coralation between any two M0 images of the different coils.
%
% This function is an lsq version of the bilinear ridge regression
% solution. this regularization can be usful to reduce the noie effect.
%
% We know that nVoxels and nCoils are redundant.  We send them in because
% we do not want to call size() every time.  We call this function a lot
% during the search.  Is there a better way to do it?
%
% AM/BW  Vistasoft Team, 2013

%% Estimate coil coeficients
G = pBasis*g;

%%
%turn your attention to the correlation (overlap) between the coil
% gain functions.  We compute the corrcoefs and and take out the lower
% triangular (non-redundant) part of these.  -1 means take everything below
% the diagonal
coefG =  tril(corrcoef(G),-1);
% Outside of this routine we already calculated
%   coefdat =tril(corrcoef(M0),-1);
% Here we compare the two correlation coefficients
% The correlation coefficients of the data (coefdat) should be larger than
% the corr coef of the gain (coefG).

corrErr = (coefdat(coefdat~=0)-coefG(coefdat~=0))./abs(coefdat(coefdat~=0)); 
% If the difference is positive, there is no regularization penalty.
corrErr(corrErr>0)=0;

%%
% Get the best PD for each position a linear sulotion
% this mske  it a nested biliner problem
% PD = zeros(nVoxels,1);
% for ii=1:nVoxels
%     PD(ii) = G(ii,:)' \ M0(ii,:)';
% end

PD=M0./G;
PD(~M0mask)=nan;
PD=nanmean(PD,2);

G  = G  .* mean(PD(PD>0));
PD = PD ./mean(PD(PD>0));

% get the predicted M0 for all of the coils
M0P = G.*repmat( PD,1,nCoils);


%The difference between measured and predicted term
SSE = (M0(M0mask) - M0P(M0mask)) ;

%SSE = (M0(:) - M0P(:))./M0(:)  ;
%%

%  the regularization term
reg=sum( W.*(g.^2));




%%
%% error
% 
%% The error is a vector with positive and negative values representing the
% M0 difference and the coil corrlations  failures.The error is PD - PDpred
%  Also the error join two other error terms one form corraltion of coils
%  and one for ridge (minimal  coefficents lenght).

err=[SSE(:)'  reg(:)'  (RegWeight)*(corrErr)'];

%err = [ M0(:) - M0P(:); (RegWeight)*(corrErr)];
end