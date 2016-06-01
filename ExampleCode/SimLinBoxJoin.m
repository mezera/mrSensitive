function [ScaleMat, LinScaleMat] = SimLinBoxJoin(Nbox, noiseLevel, Noverlap)
% Simulate data and implement linear operation to combine data from separate boxes
%
% Input:
%  Nbox -         number of local regions (boxes) 
%  noiseLevel-    noise level for the mean PD estimate in each box, a
%                 percentage of the mean.
%  Noverlap-      number of boxes that overlap a single box. The real number
%                 will be greater than the number used here. In the MRI
%                 data we typically have ~ 5000 boxes with overlap of
%                 between (0-32).  
% OutPut:
%  ScaleMat    -   Matrix of the ratio between overlapping boxes (i,j)
%  LinScaleMat -   The matrix used in the linear equation to derive the best
%                  estimate of ScaleMat given noise in the data.
%
% Notes -  Assumptions of the calculation
%
%  Estimated coil gains have been removed from each of the boxes.
%
%  The PD values of the overlapping boxes are related correctly by a single
%  scale factor. Each box has a different scale that we estimate.
%
%  We calculate the first estimate of the scale from the ratio of the data
%  in the overlapping region of each pair of boxes.
%
%  In the presence of noise, there will be some inconsistencies between
%  After the operation implemented here in which we solve a system of
%  linear equations, we derive the optimal scale factors, one for each box,
%  that provides a good global solution.
%
%  In the case of a brain, we scale all of the scale factors so that the
%  scale factor so that the PD in the ventricles is 1.0.
%
% Example:
%    nBox = 100; noiseLevel = 0.02; nOverlap = 3;
%    [S,L] = SimLinBoxJoin(nBox,noiseLevel,nOverlap);
%    mrvNewGraphWin; mesh(S);
%
% AM - Copyright Vistasoft Team, 2013
%
%% Set parameters and the unscaled box means


% Create the scalar
if notDefined('Nbox')
    Nbox    = 100;
end

if notDefined('noiseLevel')
    noiseLevel=0.05; % let say we have 5% error
end

if notDefined('Noverlap')
    Noverlap=5;
end

boxmean = round(rand(1,Nbox)*100);

%% Calculate the initial estimate of the ratio between the boxes

% Each box overlaps with few other boxes.
% We calculate the ratio of the box medians in the overlap space.
ScaleMat = zeros(Nbox,Nbox);

% In the real routine, we calculate the overlap regions so that the means
% depend only on the overlap regions
for ii=1:Nbox  % For each box
    
    % Randomly select three boxes for the overlap
    overlapBox  = randperm(Nbox);
    findBox     = find(overlapBox~=ii);
    overlapBox  = overlapBox(findBox(1:Noverlap));
    
    % Find the scalar ratios for each of the three overlap boxes with this
    % box, which is box ii
    for jj = overlapBox
        
        Ratio =( boxmean(ii)*(1+ randn(1)*noiseLevel) )./(boxmean(jj)* (1+ randn(1)*noiseLevel));
        % Ratio = boxmean(ii)./boxmean(jj); with no noise
        
        % Store the ratio of ii to jj and jj to ii
        ScaleMat(ii,jj) = 1./Ratio;
        ScaleMat(jj,ii) = Ratio;
        
    end
end

%% In actual practice
%  We only include scalars from boxes that are within some range of the hub
%  box.
%  We don't include that here.

%% Select a hub box

% We call the most connected box the hub.  If there are several with the
% same level, just take one.

% Level of connectivity by counting number of non-zero entries in each
% column of ScaleMat.
nConnect = sum(logical(ScaleMat),2);

% Find the boxes with the max number of connections.  Store the first one
% in SHub.
SHub    = find(nConnect==max(nConnect)); SHub=SHub(1);


%% Build the linear equations to solve for the all the relative scale factors

% We have removed the coil sensitivity from each of the boxes before
% getting here.  The only difference between the boxes is the unknown scale
% factor, which we created above and stored in ScaleMat.

% Therefore, the relationship between the data in the boxes is, say,
%
%   box1 = s12*box2, and perhaps box1 = s1J*boxJ
%
% We do not have an sij for every box relationship, but we do have it for
% many and an interlocking set of boxes.
%
% So, to solve for the sij given the boxes, we can set up the linear
% equation
%
%   N*box1 = \sum_j sj * boxj
%   0 = N*box1 -  \sum_j sj * boxj
%   0 = N      -  \sum_j sj * (boxj/box1)
%
% where the sum is over all the boxes that overlap with box 1.  So, we are
% going to solve for (boxj/box1).

% So, we can set up a linear equation that has a row for every box, and the
% estimated scale factors. 
LinScaleMat = zeros(Nbox,Nbox);% The key matrix for the linear relationship
LinScaleMat(SHub,SHub) = 1;    % This is the hub, with a unit scale

% For each box
for ii=1:Nbox
    % If it is not a hub
    if (ii ~=SHub)
        %
        LinScaleMat(ii,:)  = -ScaleMat(ii,:);
        LinScaleMat(ii,ii) =  nConnect(ii);
        
    end
end

%% Solve for y = [LinScaleMat] * [BoxMean]
%
% y = LinScaleMat * estBoxRatio
%
% For box 1, there is an equation 
%
%   N = 0 box11Ratio + s2 box12Ratio + s3 box13Ratio ....
%
% There is another equation for treating box 2 as the hub
%
%   N = s1 box12Ratio + 0 box22Ratio + s3 box23Ratio ....
%
% The matrix LinScaleMat contains all of the b(i,j) ratios of the median
% values. 
% We solve for the sj values, which we call the estimated box ratios.
% The estimates are exact if there is no noise.

y = zeros(Nbox,1);
y(SHub) = 1;

test = sum(LinScaleMat,1);
if ~isempty(find(test == 0, 1)) || ~isempty(find(abs(test) == Inf, 1))
    fprintf('Some boxes have no overlap. Re-run\n');
    return;
end

% Solve for the box means, assuming that the Sij are all correct.
% Given noise, the Sij can't be all correct.  So, the box means we estimate
% will not be quite right.  They will be a compromise. 
estBoxRatio = pinv(LinScaleMat'*LinScaleMat)*LinScaleMat' * y;

%% Scale value needed to correct eac box mean to agree with the hub box

% This is the correct value from the simulation
trueRatio = boxmean(SHub) ./ boxmean;

% We compare the estimated box means and the true box means
mrvNewGraphWin;
plot(estBoxRatio,trueRatio,'*')
xlabel('Estimated ratio to the hub box');
ylabel('True ratio to the hub box');
identityLine(gca); 

%% Adjust the scale factors so that the hub box has PD = 1

% MnMx=minmax(boxmean); 
% Mn=MnMx(1)*0.9;  Mx=MnMx(2)*1.1;

mrvNewGraphWin;
subplot(1,2,1);
plot(boxmean,'*')
ylabel('Original box means')
xlabel('box number')

subplot(1,2,2);
scaleRatios = estBoxRatio.*boxmean';
scaleRatios = scaleRatios / mean(scaleRatios(:));
plot(scaleRatios,'*')
ylabel('Scaled box means')
xlabel('box number')
title(sprintf('Noise level %.2f',noiseLevel));
ylim([0.8 1.2]);


%% End
