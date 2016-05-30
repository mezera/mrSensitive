%% Figure 5
%
% Illustrate howPD and coil sensativity are joined to make the M0  images %
% Code associated with Mezer, et. al. 2016, HBM
% 
%
% % AM/BW Mezer Lab & Vistaosft Team, 2013

%In this figure we plot the results of different models to estimate PD and Gain on a whole brain.
%To run the different model fit use the script saved in
% fullfile(mrSensitiveRootPath,'MethodsComparision')
% Note that since the methods recruiter additional code ( that is not our to share), here we use the saved outcome.


%%
disp('Figure 5')

%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));
%%
M0Simulated=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','M0_noise.nii.gz'));

SOS=sqrt(sum(M0Simulated.data.^2,4));

%% sub plot A
showMontage(squeeze(M0Simulated.data(:,:,45,:)));caxis([0 1500]);
colorbar off; axis off; title ('Multi coils M0')



%% sub plot B
showMontage(squeeze(SOS(:,:,45,:)));
colorbar off; axis off; title ('SOS M0')


%% sub plot c
% We show here the figures form the article with SOS coil summation and for PD free of true smooth PD in space .
% to plot to Median change to median when ever it say SOS and 
% to change to smooth PD in space change 'new' whenever we use'old'

PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','SimValues','PD.nii.gz');
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','InPut','mask.nii.gz'));
%%

PD_CoinSoS=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','COIN','GMWMmask','SOS','PDfit.nii.gz');
PD_PsedoT1SoS=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','PsedoT1','SOS','PDfit.nii.gz');
PD_unicort=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain','old','M0Files','unicort','SOS','PDfit.nii');

PD_fitmap= { PD_unicort PD_CoinSoS   PD_PsedoT1SoS};
FitName={'UNICORT' 'COIN' 'Psedo T1'};


%%
%load the fitted PD


PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;
BM=readFileNifti(BMfile); BM=logical(BM.data);

for ii=1:3
PD_fit=readFileNifti(PD_fitmap{ii});

in=PD_fit.data;
in=in.*median(PD(BM)./in(BM));
%compare to the simulation

CV = (calccod(PD(BM),in(BM))/100)
MedErr=median(abs(in(BM)-PD(BM))./PD(BM))

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
title( [' R^2= ' num2str(CV)])
grid on


Errmap=100*(in(:,:,:)-PD(:,:,:))./PD(:,:,:);
Errmap(~BM)=nan;
showMontage(Errmap(:,:,45));colormap hot
caxis([-40 40]);axis off; title(FitName(ii))
end



%%