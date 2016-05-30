function OutPut = pdBiLinearRidgeFit(M0_v,pBasis,Lambda,maxLoops,sCriterion,PD,D)
% OutPut = pdBiLinearRidgeFit_XV(M0_v,pBasis,Lambda,maxLoops,sCriterion,PD,D)
%
% This fancrtion running a while loop solving the M0=GPD bi-linear problem.
% using ridge regression with regularisation coeffisent Lamdea for the gain
% problem. im each interval of the loop. the coil gain are estimated given
% a PD values. then the PD is evaluated again given the coil gain sulotion.
% the loop stop when the solotion convarge up to the stopping criterion or
% when it reach the maxsimal iteration
%
% Input
%  M0_v:         Columns containing the 3D M0 values in a vector (nVoxels x nCoils)
%  pBasis:       Polynomial basis (nPositions x nCoef)
%  Lambda        A ridge regration regularization term
%  maxLoops      maxsimal iteration of the while loop
%  sCriterion    Stopping criterion to stop the while loop (this  define what when the problem converge)
%  PD            initial PD sulotion

% D             a wighted identity matrix that will wighting  the ridge
%                regration for each parameters D=(NcoefXnCoef)   (Ncoef=size(pBasis,2)); .
%                D is a diagonal matrix position 1,1, is wight on parameter 1 position 2,2
%                in a wight on parameter 2 ec. position (n,n) wight on parameter n.
%
% OutPut is a stracture with a list of outputs:
%  OutPut.PD             -  the final PD values
%  OutPut.Gn             -  the final coils gain values
%  OutPut.g                -  the final coils gain coeficents  (G = pBasis*g)
%  OutPut.PDchange -  a vector of the PD change in each step
%  OutPut.NumOfIter - number of iteration
%  OutPut.convarge    - 1 if convarge 0 if not
%
%   AM/BW Copyright VISTASOFT Team 2013

%% intiate parameters

if notDefined('maxLoops'), maxLoops = 1000; end
if notDefined('sCriterion') sCriterion = 1e-4;  end % Stopping criterion
   

if notDefined('PD')
    [PD]=Get_PDinit(0,[],1,M0_v,pBasis);

end


nCoils    = size(M0_v,2);
nPolyCoef = size(pBasis,2);

if notDefined('D'),
    D=eye(nPolyCoef);  D(1,1)=0; % don't work on the offset 
end

% loop and solve by ridge regration
k = 0;                % number of iteration
tryagain = 1;         % go for the while loop
PDchange = zeros(1,maxLoops);
M0change = zeros(1,maxLoops);
M0FitErr = zeros(1,maxLoops);
RelativeChange = zeros(1,maxLoops);

M0=M0_v;
%% This is the bilinear alternating solution
while tryagain == 1
    k = k+1; % count Ridge regression steps
    
    % fit linear equation coefficients using ridge regression
    g = RidgeRegressCoilfit(PD, Lambda, M0_v, pBasis,D);
    
    % Calculate the coil gain of each coil over space from the estimate
    % Calculate PD for each coil
    [PDn, Gn] = pdEstimate(M0_v, pBasis, g);
    
    
    Gn  = Gn  .* PDn(1);
    PDn = PDn ./ PDn(1);
    
    M0n = Gn.*repmat( PDn,1,nCoils);
    % Check if the new estimate differs from the one before or it's
    % converged
    PDchange(k) = std(PD - PDn);
    M0change(k) = std(M0(:) - M0n(:));
    M0FitErr(k) = std(M0_v(:)-M0n(:));
    
    % we look 100 interval back and calculate the change in fit. this will
    % allow to see if we gain any think in fitting M0.
   MimmaxErr= minmax(M0FitErr ( max(1,k-100) : k ));
   RelativeChange(k)=MimmaxErr(2)-MimmaxErr(1);
    %if PDchange(k) < sCriterion;
    if (M0FitErr(k) < sCriterion || RelativeChange(k) < sCriterion ) && k>100
        % If stable to within 1 percent, stop.
        
        % We could check the gains, rather than PD, or both
        % if std(G-Gn)<0.01
        
        % if the two solutions are stable --> stop
        tryagain=0;
    else
        % Keep going.
        
        % Update the new PD and and estimated M0
        PD = PDn;
        M0 = M0n;
        
      
            
    % We have to stop some time if it's not convarging
    if k == maxLoops
        tryagain=0;
%         if plotFlag==1
%             fprintf(' rich  % 0.f  iterations. Stop before convarge \n :',maxLoops);
%         end
    end
end



%% Outputs

% Make a list of outputs:
OutPut.PD = PDn; % the final PD values
OutPut.Gn = Gn; % the final coils gain values
OutPut.g  = g; % the final coils gain coeficents  (G = pBasis*g)
OutPut.PDchange = PDchange;  % a vector of the PD change in each step
OutPut.M0change = M0change;  % a vector of the M0 fit change in each step
OutPut.LastLoop = k;  % a vector of the PD change in each step
OutPut.M0FitErr = M0FitErr;  % a vector of the PD change in each step


OutPut.NumOfIter = k;  % number of iteration
if k == maxLoops
    OutPut.convarge=0;  % the problem while loop got to it's last interval with out convarges
else
    OutPut.convarge=1;  % the problem lo convarges
end

end

