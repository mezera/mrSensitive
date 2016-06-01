

%
%%%
%%%%%
%%%%%%
%%%%%%%%

% This is a script to the different method in figure 5 Mezer et. al. HBM 2016
% To run this script more software is needed. The software are the
% implimentation, we are testing here other approches of different authors (see below). Please ask directly
% the authors for the code when needed.

% We show here the script form the article that use SOS coil summation and for PD free of true smooth PD in space.
% to plot to Median summationchange to Median when ever it say SOS and 
% to change to smooth PD in space change 'new' whenever we use'old'
%%%%%%%
%%%%%%
%%%%
%%%
%


%% Run unicort 
% To run this section, 1. spm is needed.
% 2. unicort code is needed ( Lutti et al., MRM 2010)). Please contact Professor Nikolaus Weiskopf (UCL) 
% The output is also saved in our repository. 

 addpath(**** Path To unicort *****')

% The output is saved in our repository. See below.
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','unicort','SOS','SOS.nii');
%cd to where the result should save
cd (fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','unicort','SOS'));

%make spm stracture (spm needed)
m0=spm_vol(M0file);

%run  unicorn (aditinal code needed from the authors oforiginal article  )
[M0medGainBiasMask, P_B1] = unicort(m0,m0);

% the New PD is saved in the out put directury as hmFileName.nii

%%  Run COIN
%
%Test of COIN PD gain correction: In this file we repeat the test COIN with a full Gray and white matter mask that is used  for the ROI.

% Noterdaeme, O., Anderson, M., Gleeson, F., & Brady, S. M. (2009). Intensity correction with a pair of spoiled gradient recalled echo images. Physics in Medicine and Biology, 54, 3473?3489. doi:10.1088/0031-9155/54/11/013
% we have modified the original pipeline :
%1. We have used a white and gray matter mask to define the area that is approximately constant. The mask is segmented from the T1 map.
%2. Since COIN is 2D algorithm we have run it on the tree different plains  (axial, sagittal, coronal) and join all the 2D fits.

%We have implament this software since we could find the part of the code
%on the web.




%% Make a GM & WM mask from T1 
 T1file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','T1.nii');t1=readFileNifti(T1file);xform=t1.qto_xyz;t1=double(t1.data);
 BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));brainMask=readFileNifti(BMfile);brainMask=logical(brainMask.data);
% 
%
  GMWMmask=brainMask & t1<2.3 & t1>.6; 
  [GMWMmask] = ordfilt3D(GMWMmask,12);

  dtiWriteNiftiWrapper(single(GMWMmask),xform,fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask.nii.gz'));


%% Run COIN
 

%load a mask file for gray and white matter. 
maskfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask.nii.gz');

smoothKarnel=5;

M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','SOS.nii.gz');
outputhPath=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask','SOS')


SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel);

% Better result can be achived when the only white matter mask is used.

%% PsedoT1
%%. fit only with T1 rule volz 2012 
% Neuroimage. 2012 Oct 15;63(1):540-52. doi: 10.1016/j.neuroimage.2012.06.076. Epub 2012 Jul 14.
% Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities.
% Volz S, Nöth U, Jurcoane A, Ziemann U, Hattingen E, Deichmann R.
%Please contact Professor Ralf Deichmann , Univ. Frankfurt/Main, Germany for the code. 

addpath(**** Path to PsedoT1 Code ****)

% Note that the smooth variance PD data run was not successful 
%%
%1. get T1 with Brain mask CSF mask and M0 
T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','T1.nii.gz'));
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','mask.nii.gz'));


csffile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','csf_mask.nii.gz'));
maskfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask.nii.gz');


M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','SOS.nii');

%% set the inputs form the simulated PhantomBrain
M0=readFileNifti(M0file); xform=M0.qto_xyz;M0=double(M0.data);

T1=readFileNifti(T1file);T1=T1.data*1000;% T1 in ms
BM=readFileNifti(BMfile);BM=logical(BM.data);

GMWMmask=readFileNifti(maskfile);GMWMmask=logical(GMWMmask.data);
CSF=readFileNifti(csffile); CSF=logical(CSF.data);




%% Intiate the multi box fit parameters  1/PD=A +B/T1
A=zeros(1,8);B=zeros(1,8);
A(1) = 0.916 ; B(1) = 436; %litrature values

% a basis for estimate the new A and B.
R1basis(:,2)=1./T1(GMWMmask);
R1basis(:,1)=1;
%%

%iterate to convarge A and B 
for ii=2:8
    
myhelp=1./(A(ii-1)+B(ii-1)./T1);
% the sensativity the recive profile


% Calc. pseudoPD, allow only for pixels in maske, remove INFs and NANs, restrict to maximum of 2.0:
pseudopd=GMWMmask./(eps+myhelp);
clear myhelp
pseudopd(find(isinf(pseudopd)))=0;
pseudopd(find(isnan(pseudopd)))=0;
pseudopd=pseudopd.*(pseudopd<2);
% Calculate mask for next steps:
mymask=(pseudopd>0 & T1>100 & T1<6000).*GMWMmask;
%Now calculate quotient of M0 map and pseudopd, this should be more or less the receiver bias:
myhelp=M0./(eps+pseudopd).*mymask;
myhelp(find(isinf(myhelp)))=0;
myhelp(find(isnan(myhelp)))=0;
% This code defines an upper limit for myhelp to remove outliers:
data=myhelp(find(mymask));
mm=median(data);
ss=std(data);
for i=1:5
 data=myhelp(find(mymask & myhelp<mm+5*ss));
 mm=median(data);
 ss=std(data);
end
%Now mm+5*ss can be taken as upper limit:
myhelp=myhelp.*(myhelp<mm+5*ss);
% Create smoothed version of myhelp in all pixels of maske, this should be more or less the receiver bias field
% Details of this smoothing process to be found in above publication:
bias_field=ralf_polysmooth_3d(myhelp,GMWMmask,3,250,5,5,5,20000);
clear myhelp


% calulate PD from M0 and RP
PDi=M0./bias_field;
% scale to have PD of SCF=1
scale=mean(PDi(CSF));
PDi=PDi./scale;


% solve for A B given the new PD estiamtion
 % ( 1./PDi(GMWMmask) )= A* R1basis(:,1) + B*R1basis(:,2);

co     = R1basis \ ( 1./PDi(GMWMmask) );
A(ii)=co(1)
B(ii)=co(2)

end


PD=PDi;
PD(~BM)=nan;
%%
%the SOS case
PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','PsedoT1','PD_SOS.nii');
Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','PsedoT1','Gain_SOS.nii');



dtiWriteNiftiWrapper(single(PD),xform,PDfile);

dtiWriteNiftiWrapper(single(PD),xform,Gainfile);

