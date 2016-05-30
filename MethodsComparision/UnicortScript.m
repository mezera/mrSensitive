%% unicort test


%% Run unicort for the Median M0 image
% To run this section, 1. spm is needed.
% 2. unicort code is needed ( Lutti et al., MRM 2010)). Please contact Professor Nikolaus Weiskopf for the code(UCL) 

% addpath(Path To unicort')

% The output is saved in our repository. See below.
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','unicort','sos','SOS.nii');
%cd to where the result should save
cd (fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','unicort','sos'));

%make spm stracture (spm needed)
m0=spm_vol(M0file);

%run  unicorn (aditinal code needed from the authors oforiginal article  )
[M0medGainBiasMask, P_B1] = unicort(m0,m0);

%% Run unicort for the SOS M0 image
% To run this section, 1. spm is needed.
% 2. unicort code is needed ( Lutti et al., MRM 2010)). Please contact Professor Nikolaus Weiskopf (UCL) 
% The output is saved in our repository. See below.

M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','unicort','sos','sos.nii');
%cd to where the result should save
cd (fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','unicort','sos'));

%make spm stracture
m0=spm_vol(M0file);

%run  unicorn
[M0medGainBiasMask, P_B1] = unicort(m0,m0);

%% smooth PD sos
% to try the M0 with some true smooth PD on it evaluate this line instade
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','unicort','sos','sos.nii');
%cd to where the result should save
cd (fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','unicort','sos'));


%make spm stracture
m0=spm_vol(M0file);
%run  unicorn
[M0sosGainBiasMask, P_B1] = unicort(m0,m0);

%% smooth PD median
% to try the M0 with some true smmmoth PD on it evaluate this line instade
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','unicort','median','M0_median.nii');
%cd to where the result should save
cd (fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','unicort','median'));

%make spm stracture
m0=spm_vol(M0file);

%run  unicorn
[M0medGainBiasMask, P_B1] = unicort(m0,m0);
