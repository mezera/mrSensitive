function err = errFitRidgeNestBiLinear_lsq_v1(g,M0,pBasis,nCoils,W,M0mask)
%Error function for bilinear coil gain estimate for fminsearch
%
%  err = errFitRidgeNestBiLinear(g,M0,pBasis,nVoxels,nCoils,W,FitMask)
%
% g       - Coil gain polynomial coefficients (nCoeffs x nCoils)
% M0      - Measurement (nVoxels x nCoils)
% pBasis  - Polynomial basis in columns (nVoxels x nCoeffs)
% nVoxels - Number of voxels
% nCoils  - Number of measurement coils
% W       - The ridge regression coefficients
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

% Get the best PD for each position a linear sulotion
% this mske  it a nested biliner problem
PD=M0./G;
PD(~M0mask)=nan;
PD=nanmean(PD,2);



G  = G  .* mean(PD(PD>0));
PD = PD ./mean(PD(PD>0));

% get the predicted M0 for all of the coils
M0P = G.*repmat( PD,1,nCoils);




%  the regularization term
reg=sum( W.*(g.^2));
%reg= W.*g;

%reg=( W.*(g.^2));
%The difference between measured and predicted term
SSE = (M0(M0mask) - M0P(M0mask)) ;
%SSE = (M0(:) - M0P(:))./M0(:)  ;

% join the two error terms
err=[SSE(:)'  reg(:)'];

end


% PD1 = zeros(nVoxels,1);
% for ii=1:nVoxels
%     mask=M0mask(ii,:);
%     PD1(ii) = G(ii,mask)' \ M0(ii,mask)';
% end

