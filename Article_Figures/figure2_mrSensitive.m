%% figure 2
clear
disp('Figure 2')

% Making a simulated data by PD and homogeneous phantom data M0.
% Code associated with Mezer, et. al. 2016, HBM
%
% % AM/BW  Mezer Lab Wandell Lab & Vistaosft Team, 2013


% Also see Sim_M0_Noise.m in: fullfile(mrSensitiveRootPath,'MethodsComparision')
%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));
%%
PD=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','PD.nii.gz'));

M0CoilSensitivity=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','Gain.nii.gz'));
% this coil image masked by the brain voxels
%The raw coil M0 is here:
%fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0'));

M0Simulated=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','M0_noise.nii.gz'));

PD_T1reg=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','T1reg','PDfit.nii.gz');

PDestimate=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','T1reg','PDfit.nii.gz'));

%%
showMontage(PD.data(:,:,45,1));
colorbar off; axis off; title ('PD')


showMontage(M0CoilSensitivity.data(:,:,45,1));
colorbar off; axis off; title ('Sensitivity')

showMontage(M0Simulated.data(:,:,45,1));
colorbar off; axis off; title ('M0')

showMontage(PD.data(:,:,45,1));
colorbar off; axis off; title ('Estimated PD')







%% 
% Noise=M0Simulated.data(:,:,:,1)./PD.data-M0CoilSensitivity.data(:,:,:,1);
% Mask=PD.data>0;
% Noise(~Mask)=0;
% showMontage(Noise(:,:,45));
% colorbar off; axis of; title ('Noise')


