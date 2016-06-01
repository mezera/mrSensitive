%% figure 7
% 
% Code associated with Mezer, et. al. 2016, HBM
% 
%
% % AM/BW Mezer Lab & Vistaosft Team, 2015
%%

clear
disp('Figure 7')

% Compares the same brain measured using 8 and 32 channel coil data
%
% To test the reliability of our method we scan the same subject in two
% different scans with 8 and 32 channels coils. we compare the PD estimate
% in the two scans

% Raw Data can be found here as nifti:  http://purl.stanford.edu/nn554zr6949.
% In there scan-1 is the 8ch and scan-2 is the 32ch.

% AM/BW Vistaosft Team, 2013

% Running the solution  for PD from raw image will take many hours.
% For the analysis of PsedoT1 and Uincort see details in figure5_mrSensitive.m.
%Note that since these methods recruiter additional code ( of other authors), here we use the saved outcome.
% For running the correlation method use  https://github.com/mezera/mrQ VS.1
% For running the Local-psedoT1 approach use https://github.com/mezera/mrQ VS. 2

%% Loading the calculated results



%% masking the area for comparison. 
%Brain masks:

BMS1=(fullfile(mrSensitiveRootPath,'ExampleData','Ch8Ch32SavedMaps' ,'Masking' ,'8','brainMask.nii.gz'));
BMS2=(fullfile(mrSensitiveRootPath,'ExampleData','Ch8Ch32SavedMaps' ,'Masking' ,'32','brainMask.nii.gz'));

BMS1=readFileNifti(BMS1);
BMS2=readFileNifti(BMS2);

% load to make a GM WM mask using T1 values

T1S1=(fullfile(mrSensitiveRootPath,'ExampleData','Ch8Ch32SavedMaps' ,'Masking' ,'8','T1.nii.gz'));
T1S2=(fullfile(mrSensitiveRootPath,'ExampleData','Ch8Ch32SavedMaps' ,'Masking' ,'32','T1.nii.gz'));

T1S1=readFileNifti(T1S1);
T1S2=readFileNifti(T1S2);

%% make a brain mask that containing only brain tissue 

BM = logical(BMS2.data) & logical(BMS1.data);

BM(:,:,1:50)   = 0;
BM(:,:,140:end)= 0;

% To assure we are in GM or WM we want T1 that is greater then 0.5 and smaller then 2.5. % 
BM1=BM & T1S1.data>0.6 & T1S1.data<2.5 & T1S2.data>0.6 & T1S2.data<2.5 ;


%% PD output maps

Maps_path= fullfile(mrSensitiveRootPath,'ExampleData','Ch8Ch32SavedMaps' );


MethodName={'Unicort','Corr','PsedoT1','T1Reg','LocalPsedoT1'};

%%
close all

clear CV MAPE
%%
for ii=1:5
% PD files:

WF1=readFileNifti(fullfile(Maps_path,'32',MethodName{ii},'PDfit.nii.gz'));
WF2=readFileNifti(fullfile(Maps_path,'8',MethodName{ii},'PDfit.nii.gz'));


%%


%% let's remove any constant different between the PD maps
WF2.data=WF2.data*median(WF1.data(BM1)./WF2.data(BM1));


%% The coefficient of determination (R^2)
CV(ii) = (calccod(WF1.data(BM1),WF2.data(BM1))/100);


%% Panel 7B Plot a  2D histogram of the error
 

nBins = 155;
[n,x,y] = mrQ_hist2d(WF1.data(BM1),WF2.data(BM1),nBins);

maxN = ceil(max(n(:))/10)*10;

mrvNewGraphWin;
image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
identityLine(gca);
axis square xy; axis image;
xlim([0.6 1]); ylim([0.6 1])
ylabel('PD  32ch' ,'FontSize',16);
xlabel('PD  8ch' ,'FontSize',16);
set(gca,'FontSize',16)
title( [ MethodName{ii} ' R^2= ' num2str(CV(ii)) ])
grid on

if ii==5; text(0.62,0.92,'Not shown in the article ');end


%% Panel 7c Error maps  

Err=(WF1.data-WF2.data)./((WF2.data+WF2.data)./2);
Err(isnan(Err))=0;Err(isinf(Err))=0;
Err(~BM)=-100;
mrvNewGraphWin;
imagesc(Err(:,:,90)*100); colormap hot; caxis([-10 10]); axis off
MAPE(ii)=mean(abs(Err(BM))) ;
title([MethodName{ii} ' Precent error= ' num2str(MAPE(ii)) ])
 
if ii==5; text(10,200,'Not shown in the article ');end

end


%% Panel 7A PD maps
showMontage(WF1.data(:,:,90,1));
colorbar off; axis off; title ('PD 8ch'); caxis([0 1.2])
showMontage(WF2.data(:,:,90,1));
colorbar off; axis off; title ('PD 32ch');caxis([0 1.2])
