%% figure 2
clear
disp('Figure 2')

% Making a simulated data by multiple PD and homogeneous phantom data M0.
% Code associated with Mezer, et. al. 2016, HBM
%
% % AM/BW  Mezer Lab Wandell Lab & Vistaosft Team, 2013


% Also see Sim_M0_Noise.m in: fullfile(mrSensitiveRootPath,'MethodsComparision')
%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));
%%
PD=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','PD.nii.gz'));

M0CoilSensativity=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','Gain.nii.gz'));
% this coil image masked by the brain voxels
%The raw coil M0 is here:
%fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0'));

M0Simulated=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','M0_noise.nii.gz'));

%%
showMontage(PD.data(:,:,45,1));
colorbar off; axis off; title ('PD')


showMontage(M0CoilSensativity.data(:,:,45,1));
colorbar off; axis off; title ('Sensitivity')

showMontage(M0Simulated.data(:,:,45,1));
colorbar off; axis off; title ('M0')









%% 
% Noise=M0Simulated.data(:,:,:,1)./PD.data-M0CoilSensativity.data(:,:,:,1);
% Mask=PD.data>0;
% Noise(~Mask)=0;
% showMontage(Noise(:,:,45));
% colorbar off; axis of; title ('Noise')


