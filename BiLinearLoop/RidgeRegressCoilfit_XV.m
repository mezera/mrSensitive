function g = RidgeRegressCoilfit_XV(PD, Lambda, M0, pBasis,D,Mask,W)
% Solve for gain parameters using ridge regression
% 
%  g = RidgeRegressCoilfit(PD,Lamda,M0,Pbasis)
%
% PD:                 Vector of estimated proton densities at nPositions
% Lambda:        Weight for the ridge
% M0:                The coil data, nPositions x nCoils
% pBasis:           Polynomial basis
% D                    a wighted identity matrix that will wighting  the ridge
%                        regration for each parameters D=(NcoefXnCoef)   (Ncoef=size(pBasis,2)); .
%                        D is a diagonal matrix position 1,1, is wight on parameter 1 position 2,2
%                        i n a wight on parameter 2 ec. position (n,n) wight on
%                        parameter n.  Defult all wight are equal 1.
%W                     A diagonal matrix  (nVoxels X nVoxels)  of the wight
%                         of each raw im M0 for the lsq fit
%
% AM/BW VISTASOFT 2013


% please let's verify that this is right
Phat = diag(PD)*pBasis;

if notDefined('D')
    D=eye(size(Phat,2));
end
if notDefined('W')
    wightedFlag=false;
else
    wightedFlag=true;
end
%g1 = (Phat'*Phat + Lambda*D)^-1 * Phat'*M0;

for ii=1:size(Mask,2)
CMask=Mask(:,ii);
g(:,ii) = (Phat(CMask,:)'*Phat(CMask,:) + Lambda*D)^-1 * Phat(CMask,:)'*M0(CMask,ii);

end



if wightedFlag
    %wighted lsq
    %wighting by signal higer signal have higer SNR as noise is similar in space 
    for ii=1:size(Mask,2)
        CMask=Mask(:,ii);
        g(:,ii) = (Phat(CMask,:)'*W*Phat(CMask,:) + Lambda*D)^-1 * Phat(CMask,:)'*W*M0(CMask,ii);
    end
end
