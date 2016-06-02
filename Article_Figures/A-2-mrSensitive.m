%Figure A-2
clear
disp('Figure A-2')

%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));


%% run the simulation many times. with differnt conditions (noise,coil sensativity, pd 3D stracture)

% we will run the simulation many conditions (5 differnt sensativties function; 10 PD 3D stracture; 10 reapets of random noise =1250 ).
%we will solve each case problem with multipal fiiting approches. ALS(alternating lease squre, ridgre regression, corralation regularizaton, T1 regularization) 
%This take long time !!
fprintf('These simulations and fits are long (~30 hours). We will skip it and load the precalculate results. \n To run the simulations please evaluate lines (22-44). \n')

% you may load the saved output below in the section Bellow

% Pleae evaluate this:
%
% %the order of the pholinomyal that are used to simulate and fit the coil sensativity.
% % When CSorder=0 simulate order and fit order are the same. when CSorder is 1 the fit order will be order+1 from  the simulate order. ect.
%         CSorder=0;%
% 
% 
%         noise=2; %noise level
% 
%         loci=[5 4 3 2 1]; % phantom coil locations
%  % repeat ii time each with different random noise 
%  
%          ModelTypes=[1 0 1 1 1 0]; 
%   % we can solve the problem with differnt approch [Correlation regularization,[],ridge regresion,T1 regularization,ALS,multi T1 tissue type regularization]
%  
%         for ii=1:25
%             ii
%             % run ovr different phantom position and use it to simulate the coil gains
%             for loc=loci 
%                 % run over different PD stracture in space
%                 for PD=1:10 
%                     PD
%               [ Err(ii,loc,PD,:)]=A-2Simulations(noise,loc,num2str(PD-1),ModelTypes,CSorder) ;     
%               save( fullfile(mrSensitiveRootPath,['MultiRegFitsOnOrderNoise_' date]),'Err')
%                 end
%             end
%         end


%% The precalculate results can be load here

%to load evaluate this line:

load (fullfile(mrSensitiveRootPath,'ExampleData','simulationResults' ,'BLRegFitsONOrderNoise_01-Mar-2014.mat'))

%% The bar graph figures

mrvNewGraphWin('Figure A-2','wide');

%% Figure 4A
% the figure  use the Err output of figure4Simulations.m calculations.
%For each simulation and fits we save theresidual between the estimate and simulte PD.
%For each simulation and fits we look on the residual between the estimate and simulte PD.
%The error is  de-meaned becouse ther eis alqay a constat free pparameter that we can't fit



%   Each simulation is fited with different approch. each fitting method mean absulot error is saved: 

% 1. ALS (no regularization approch
N0=[Err(:,:,:).Err_Noreg];

% 2. corralation regularization approch
Cor=[Err(:,:,:).Err_cor];

%3. the ridge regression approch
Rig=[Err(:,:,:).Err_ridgeBL_XV];

% 4. the T1 regularization approch
T1=[Err(:,:,:).Err_T1reg];

% we will plot a box plot of the Log 10 of the mean absulot error.
subplot(1,2,1) 
boxplot(log10([ N0(2:2:end)'  Cor(2:2:end)' Rig(2:2:end)' T1(2:2:end)' ]),...
     'notch','on', 'labels',{'Non Reg', 'Cor', 'Ridge',   'T1'},'datalim',[0 3.6],'symbol', 'w.' ,'outliersize',1e-8,'colors','k');
 ylabel('Log10 ( mean percent error ) ','FontSize',16)
set(gca,'FontSize',17,'ylim',[0 3.5]); 
title('Figure A-2A','FontSize',16)
set(gca,'YTick',[0.5 1 1.5 2 2.5 3])
set(gca,'YTickLabel',{'0.5' '1' '1.5' '2' '2.5' '3'})
set(gca,'XTick',[ 1 2 3 4])
set(gca,'XTickLabel',{'None' 'Correlation' 'Tikonov' 'T1'})
grid on
%% Figure A-2B The Spatial coherence of the errors
% the figure  use the Err output of figure4Simulations.m calculations.
%For each simulation and fits we save theresidual between the estimate and simulte PD.
%For each simulation and fits we look on the residual between the estimate and simulte PD.
%The error is  de-meaned becouse ther eis alqay a constat free pparameter that we can't fit
%In order to test for the exsistence of spatial bias (residual spatial coherence), we have also fit a linear model to the (error) residual in space. 
%We will use the (r) pirsion correlation coeecisentof the lineal model to the
%residual . the variance explain (r^2) by the linear model of the error is
%the is a measure of the spatial bias (residual spatial coherence).


%Each simulation is fited with different approch. each fitting method residual spatial coherence is saved (we saved r so we squre it now): 

% 1. ALS (no regularization approch
N0_ve=[Err(:,:,:).Err_Noreg_ErrResid].^2;

% 2. corralation regularization approch
Cor_ve=[Err(:,:,:).Err_cor_LinErrResid].^2;

%3. the ridge regression approch
Rig_ve=[Err(:,:,:). Err_ridgeBL_LinErrResid].^2;

% 4. the T1 regularization approch
T1_ve=[Err(:,:,:).Err_T1reg_ErrResid].^2;

% we will plot a box plot of the residual spatial coherence 
subplot(1,2,2) 
boxplot(([ N0_ve(:) Cor_ve(:) Rig_ve(:)  T1_ve(:) ]),...
     'notch','on', 'labels',{'No Reg' ,'Cor' ,'Ridge',  'T1'}, 'symbol', 'w.' ,'outliersize',1e-8,'colors','k');
 ylabel(' Spatial coherence of the errors ','FontSize',16)
title('A-2B','FontSize',16)
set(gca,'FontSize',17); 
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1 ])
set(gca,'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1'})
set(gca,'XTick',[ 1 2 3 4])
set(gca,'XTickLabel',{'None' 'Correlation' 'Tikonov' 'T1'})

grid on

