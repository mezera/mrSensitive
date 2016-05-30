%% The reason for the work and testing figure
%%
disp('Test M0 methods using T1 global')


%%. fit only with T1 rule volz 2012 
% Neuroimage. 2012 Oct 15;63(1):540-52. doi: 10.1016/j.neuroimage.2012.06.076. Epub 2012 Jul 14.
% Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities.
% Volz S, Nöth U, Jurcoane A, Ziemann U, Hattingen E, Deichmann R.
%Please contact Professor Ralf Deichmann , Univ. Frankfurt/Main, Germany for the code. 

addpath('/Users/avivmezer/Documents/articles/PD_methodArticle/HBM/Weiskopf_Deichmann_Code/Ralf_T1PD')


%%


%1. get T1 with Brain mask CSF mask and M0 
T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','T1.nii.gz'));
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','mask.nii.gz'));


csffile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','csf_mask.nii.gz'));
maskfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask.nii.gz');

%WE can try Median and SOS coil summing. and PD with and without smoothbias 
%M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','M0_median.nii');
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','SOS.nii');
% the smoothbiased version
%M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0Files','SOS.nii');
%M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0Files','M0_median.nii');

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
%median
% PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','PD_median.nii');
% Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','Gain_median.nii');

%the SOS case
PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','PD_SOS.nii');
Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','Gain_SOS.nii');

% smooth
% PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','PsedoT1','PD_median.nii');
% Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','PsedoT1','Gain_median.nii');
% 
% %the SOS case
% PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','PsedoT1','PD_SOS.nii');
% Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','PsedoT1','Gain_SOS.nii');


dtiWriteNiftiWrapper(single(PD),xform,PDfile);

dtiWriteNiftiWrapper(single(PD),xform,Gainfile);



%%


Gfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','Gain.nii.gz');

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

% median
% PDpsedofile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','PD_median.nii');
% OutPutFileSos=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','Medoutputs');

%sos
PDpsedofile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','PD_SOS.nii');
OutPutFileSos=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','PsedoT1','outputs');

%%
PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;
PD_t1=fullfile(PDpsedofile);in=readFileNifti(PD_t1);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
MedErr=median(abs(in(BM)-PD(BM))./PD(BM))
save(OutPutFileSos,'CV','MedErr')

nBins = 155;
[n,x,y] = mrQ_hist2d(PD(BM),in(BM),nBins);
maxN = ceil(max(n(:))/10)*10;

mrvNewGraphWin;
image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
identityLine(gca);
axis square xy; axis image;
xlim([0.5 1.1]); ylim([0.5 1.1])
ylabel('PD  sim' ,'FontSize',16);
xlabel('PD  fit' ,'FontSize',16);
set(gca,'FontSize',16)
title( [' R^2= ' num2str(CV)])
grid on



Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off

