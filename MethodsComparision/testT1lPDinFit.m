%% The reason for the work and testing figure
%%
disp('Test M0 methods ')


%%. fit only with T1 rule volz 2012 
% Neuroimage. 2012 Oct 15;63(1):540-52. doi: 10.1016/j.neuroimage.2012.06.076. Epub 2012 Jul 14.
% Quantitative proton density mapping: correcting the receiver sensitivity bias via pseudo proton densities.
% Volz S, Nöth U, Jurcoane A, Ziemann U, Hattingen E, Deichmann R.



% %% simulate with PD fitted with Corr
% % make a mrQ directory and parameters
% M0file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'new' ,'InPut','M0_noise.nii.gz'));
% %1. make M0 with G PD and noise have T1
% T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'new','InPut','T1.nii.gz'));
% BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'new', 'InPut','mask.nii.gz'));

%dirpath ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting'


%% %% Get the data ( the data was simulate from the fitted T1 and fitted PD with T1reg)

%        make a M0 from multi coils

%M0file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old' ,'InPut','M0_noise.nii.gz'));
% M0=readFileNifti(M0file);
% M0=mean(M0.data,4); % the mean is good aprocsimate of the PD so we will start with that.
% save the mean M0
% M0file=fullfile(pwd,'M0');
% dtiWriteNiftiWrapper(single(M0),xform,M0file);
% M0file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'new' ,'InPut','M0_noise.nii.gz'));
%  M0=readFileNifti(M0file);
%  M0=mean(M0.data,4); % the mean is good aprocsimate of the PD so we will start with that.
% % save the mean M0
%  M0file=fullfile(pwd,'M0n');
%  dtiWriteNiftiWrapper(single(M0),xform,M0file);




%1. get T1 with Brain mask CSF mask and M0 
T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','T1.nii.gz'));
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','mask.nii.gz'));


csffile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','csf_mask.nii.gz'));

dirpath ='/home/aviv.mezer/Documents/Reaserch/PD_article/SimResults'
%dirpath ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting'
M0file=fullfile(dirpath,'M0');
%% set the inputs form the simulated PhantomBrain
M0=readFileNifti(M0file);

xform=M0.qto_xyz;
M0=double(M0.data);

T1=readFileNifti(T1file);T1=T1.data*1000;% T1 in ms
BM=readFileNifti(BMfile);BM=logical(BM.data);
GMWMmask=BM & T1<2300 & T1>600; 

CSF=readFileNifti(csffile); CSF=logical(CSF.data);

%% Intiate the multi box fit parameters  1/PD=A +B/T1
A=zeros(1,7);B=zeros(1,7);
A(1) = 0.916 ; B(1) = 436; %litrature values

% a basis for estimate the new A and B.
R1basis(:,2)=1./T1(GMWMmask);
R1basis(:,1)=1;

%iterate to convarge A and B 
for ii=2:7
    
PDp=1./(A(ii-1)+B(ii-1)./T1);
% the sensativity the recive profile
RPp=M0./PDp;

% smooth the estimation of GM & WM
[RPi]=mrQ_smoothGain(RPp,GMWMmask);

% calulate PD from M0 and RP
PDi=M0./RPi;
% scale to have PD of SCF=1
scale=mean(PDi(CSF));
PDi=PDi./scale;


% solve for A B given the new PD estiamtion
 % ( 1./PDi(GMWMmask) )= A* R1basis(:,1) + B*R1basis(:,2);

co     = R1basis \ ( 1./PDi(GMWMmask) );
A(ii)=co(1)
B(ii)=co(2)

end

PDp=1./(A(ii)+B(ii)./T1);
% the sensativity the recive profile
RPp=M0./PDp;

% smooth the estimation of GM & WM
[RPi]=mrQ_smoothGain(RPp,GMWMmask);

% calulate PD from M0 and RP
PD=M0./RPi;
% scale to have PD of SCF=1
scale=mean(PDi(CSF));
PD=PD./scale;

PD(~BM)=nan;

PDfile=fullfile(dirpath,'PD');
dtiWriteNiftiWrapper(single(PD),xform,PDfile);
Gainfile=fullfile(dirpath,'Gain');
dtiWriteNiftiWrapper(single(PD),xform,Gainfile);


%% Run with the code of Ralf Deichmann, Univ. Frankfurt/Main, Germany, May 2014
% this code dir hold private code of Ralf Deichmann
cd /home/avivm/Documents/PD_article/R1_Ralf

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

% pd_map=100*m0./(eps+abs(bias_field)).*GMWMmask;
% pd_map(find(isinf(pd_map)))=0;
% pd_map(find(isnan(pd_map)))=0;
% % PD should not exceed 100, but this is only an approximate PD, so chose limit of 200:
% pd_map=pd_map.*(pd_map>=0 & pd_map<200)+200*(pd_map>=200);
% 
% 

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

PDfile=fullfile(dirpath,'PD_ralf');
dtiWriteNiftiWrapper(single(PD),xform,PDfile);
Gainfile=fullfile(dirpath,'Gain_ralf');
dtiWriteNiftiWrapper(single(PD),xform,Gainfile);



%%


Gfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','Gain.nii.gz');

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_t1=fullfile(dirpath,'PD');

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;
%%
in=readFileNifti(PD_t1);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot

%%
PD_ralf=fullfile(dirpath,'PD_ralf');

in=readFileNifti(PD_ralf);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot

%%


PD_t1local='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/PDFItMethodTesting_multiCondisions/mrQ_R_5_C_1_I_5/PD.nii.gz';

in=readFileNifti(PD_t1local);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot


%%
Gfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','Gain.nii.gz');

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_t1=fullfile(dirpath,'PD');

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;
%%
in=readFileNifti(PD_t1);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot

%%
PD_ralf=fullfile(dirpath,'PD_ralf');

in=readFileNifti(PD_ralf);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot
%% local fit
PD_t1='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/PDFItMethodTesting_multiCondisions/mrQ_R_5_C_1_I_5/PD.nii.gz';

%%
in=readFileNifti(PD_t1);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot

%% new
%%
Gfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','SimValues','Gain.nii.gz');

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_t1='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/mrQ_R_5_C_1_I_6/PD.nii.gz';

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;
%%
in=readFileNifti(PD_t1);

in=in.data;
in=in.*median(PD(BM)./in(BM));

CV = (calccod(PD(BM),in(BM))/100)
median(abs(in(BM)-PD(BM))./PD(BM))

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
showMontage(Errmap(:,:,45));

caxis([-50 50]); colormap hot