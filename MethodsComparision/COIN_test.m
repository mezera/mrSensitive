%%
%Test of COIN PD gain correction: In this file we repeat the test COIN with a full Gray and white matter mask that is used  for the ROI.

% Noterdaeme, O., Anderson, M., Gleeson, F., & Brady, S. M. (2009). Intensity correction with a pair of spoiled gradient recalled echo images. Physics in Medicine and Biology, 54, 3473?3489. doi:10.1088/0031-9155/54/11/013
% we have modified the original pipeline :
%1. We have used a white matter mask to define the area that is approximately constant. The mask is segmented from the T1 map.
%2. Since COIN is 2D algorithm we have run it on the tree different plains  (axial, sagittal, coronal) and join all the 2D fits.


%% use the saved WM mask or recalculate it (bellow) the segmented T1 map or evaluate the code below



%% Make a WM mask from T1
% T1file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','T1.nii');t1=readFileNifti(T1file);xform=t1.qto_xyz;t1=double(t1.data);
% BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));brainMask=readFileNifti(BMfile);brainMask=logical(brainMask.data);
% 
%
%  GMWMmask=brainMask & t1<2.3 & t1>.6; 
%  [GMWMmask] = ordfilt3D(GMWMmask,12);
%  dtiWriteNiftiWrapper(single(GMWMmask),xform,fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask.nii.gz'));


%% Run COIN
%calculate coil gains from the PD_fit and M0. M0 is 4D image the 4th
%dimention is the image of each coil
% M0=Gain X PD  --> Gain= M0 / PD;
% PD fit is not full image. we have values only where the boxes fit where done
% and got a "good sulotion". hopfully we have most of the image.
% It is given that the gain function vart smoothly in space.
% Therefore to fill the area where PD data is missing we will model the
% Gain as a smooth function in space.

%% you may run  the line below and run COIN  or load the the saved maps (see  below)


maskfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask.nii.gz');

smoothKarnel=5;



%%

% COIN input SOS M0 file
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','SOS.nii.gz');
outputhPath=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','SOS')
%run COIN on WM mask (the idea that any smooth variation in M0 over the mask is coil bias )


SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel);


%COIN  input Median M0 file
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','M0_median.nii.gz');
outputhPath=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','Median')
%run COIN on WM mask (the idea that any smooth variation in M0 over the mask is coil bias )
SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel);


%% to try the M0 with some true smooth PD on it evaluate this line instade
% COIN input SOS M0 file
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','SOS.nii.gz');
outputhPath=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','SOS');
%run COIN on WM mask (the idea that any smooth variation in M0 over the mask is coil bias )
SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel);


%COIN  input Median M0 file
M0file=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','M0_median.nii.gz');
outputhPath=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','Median');
%run COIN on WM mask (the idea that any smooth variation in M0 over the mask is coil bias )
SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel);


       
        
%%
PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_CoinSoS=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','SOS','PDfit.nii.gz');
PD_CoinMed=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','Median','PDfit.nii.gz');

OutPutFileSos=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','SOS','outputs');
OutPutFileMed=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0files','COIN','GMWMmask','Median','outputs');



%% make figures for SOS M0
%load the fitted PD

PD_Coin=readFileNifti(PD_CoinSoS);

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

in=PD_Coin.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation

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


%calculate the bias map

Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off
%% make figures for Medinan M0
%load the fitted PD
PD_Coin=readFileNifti(PD_CoinMed);

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

in=PD_Coin.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation
CV = (calccod(PD(BM),in(BM))/100)
MedErr=median(abs(in(BM)-PD(BM))./PD(BM))
save(OutPutFileMed,'CV','MedErr')

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

%calculate the bias map

Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45)); colormap hot
caxis([-40 40]);axis off

%% the smooth case
PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','InPut','mask.nii.gz'));
BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_CoinSoS=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','SOS','PDfit.nii.gz');
PD_CoinMed=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','Median','PDfit.nii.gz');
OutPutFileSos=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','SOS','outputs');
OutPutFileMed=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','new','M0files','COIN','GMWMmask','Median','outputs');
%% make figures for SOS M0
%load the fitted PD

PD_Coin=readFileNifti(PD_CoinSoS);

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

in=PD_Coin.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation

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


%calculate the bias map

Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off
%% make figures for Medinan M0
%load the fitted PD
PD_Coin=readFileNifti(PD_CoinMed);

PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

in=PD_Coin.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation
CV = (calccod(PD(BM),in(BM))/100)
MedErr=median(abs(in(BM)-PD(BM))./PD(BM))
save(OutPutFileMed,'CV','MedErr')

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

%calculate the bias map

Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45)); colormap hot
caxis([-40 40]);axis off
