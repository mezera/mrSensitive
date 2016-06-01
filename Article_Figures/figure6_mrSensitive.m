%% Figure 6
%
% 
%
% % AM/BW  Mezer Lab Wandell Lab & Vistaosft Team, 2014

%In this figure we plot the results of different models to estimate PD and Gain on a whole brain.
%To run the different model fit use the script saved in
% fullfile(mrSensitiveRootPath,'MethodsComparision'). MultiMethodScript_figure6.m
% Note that since the fits take very long. here we use the saved outcome.


%%
disp('Figure 6')

%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));


% We show here the figures form the article with PD that is simulated free of smooth PD values in space .
% to change to smooth PD in space change 'new' whenever we use 'old'

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
%%


PD_T1reg=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','T1reg','PDfit.nii.gz');
PD_Tikonov=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','Corr','PDfit.nii.gz');
PD_corr=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','Ridge','PDfit.nii');

PD_fitmap= {  PD_corr   PD_Tikonov PD_T1reg};
FitName={'Tikonov','Corr' 'T1reg' };


%%
%load the fitted PD


PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

for ii=1:3
    BM=readFileNifti(BMfile); BM=logical(BM.data);

PD_fit=readFileNifti(PD_fitmap{ii});

in=PD_fit.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation

CV = (calccod(PD(BM),in(BM))/100);
MedErr=median(abs(in(BM)-PD(BM))./PD(BM));

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
xlabel(['PD  fit ' FitName(ii) ],'FontSize',16);
set(gca,'FontSize',16)
grid on
title( [' R^2= ' num2str(CV)])


Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off; title([ FitName(ii) 'NAPE ' num2str(MedErr)])
end



%% Not a figure %%
% to plot to Median change to 'Median' where is say 'SOS'  
% to plot to the smooth PD change to 'new' here is say 'old'

PD_fit=readFileNifti(fullfile(mrSensitiveRootPath,'/ExampleData/PhantomBrain/old/M0Files/LocalT1/SOS/PDfit.nii.gz'));
    BM=readFileNifti(BMfile); BM=logical(BM.data);




in=PD_fit.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation

CV = (calccod(PD(BM),in(BM))/100);
MedErr=median(abs(in(BM)-PD(BM))./PD(BM));

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
xlabel(['PD  fit local-T1 -  not shown in the article' ],'FontSize',16);
set(gca,'FontSize',16)
grid on
title( [FitName(ii) ' R^2= ' num2str(CV) ])


Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off; title([ FitName(ii) 'NAPE ' num2str(MedErr) ' not shown in the article'] )

