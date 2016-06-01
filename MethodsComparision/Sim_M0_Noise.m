%1. make M0 with G PD and noise have T1
T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'InPut','T1.nii.gz'));
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'InPut','mask.nii.gz'));



PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','PD.nii.gz');
Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','Gain.nii.gz');
G=readFileNifti(Gainfile);G=G.data;
PD=readFileNifti(PDfile);xform=PD.qto_xyz;PD=PD.data;

nVoxels=size(PD);
noiseLevel=2;

for ii=1:size(G,4)
M0(:,:,:,ii)=PD.*G(:,:,:,ii) +  randn(nVoxels)*noiseLevel;
end
M0file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'InPut','M0_noise.nii.gz'));

dtiWriteNiftiWrapper(single(M0),xform,M0file);



%%
% % load the example data from http://purl.stanford.edu/qh816pc3429 and use
% % the PD image
% T=readFileNifti('/home/avivm/Documents/PD_article/PD_example/WF_map.nii.gz');
% % let's resample to 2 2 2
% outMm=[ 2 2 2];
% sz=size(T.data);
% bb = mrAnatXformCoords(T.qto_xyz, [1 1 1; sz(1:3)]);
% interp=1;
% [PDresamp, M0UnderSamp_Xform] = mrAnatResliceSpm(double(T.data), inv(T.qto_xyz), bb, outMm, interp, 0);
% 
% %% we will make an M0 by multiple the PD and the M0 of the phantom.
% %make head mask and mask 
% load /biac4/wandell/biac2/wandell2/data/WMDevo/test/20111020_1294/SPGR_1/Align_0.9375_0.9375_1/dat_aligned.mat
% [Dat, M0UnderSamp_Xform] = mrAnatResliceSpm(double(s(1).imData), inv(T.qto_xyz), bb, outMm, interp, 0);
% 
% BM=readFileNifti('/biac4/wandell/biac2/wandell2/data/WMDevo/test/20111020_1294/SPGR_1/Align_0.9375_0.9375_1/BMresmp_2_2_2.nii.gz');
% mask=BM.data>0 & PDresamp>0;
% 
% %bring coil gain from phantom
% 
% 
% %% make a M0 data
% 
% 
% Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'SimValues','Gain.nii.gz');
% G=readFileNifti(Gainfile);G=G.data;
% 
% 
% Gsos=sqrt(sum(G.^2,4));
% 
% 
% 
% 
% 
% M0=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0','AllCoilsM0_phantomExample.nii.gz'));
% M0=sqrt(sum(M0.data.^2,4)); % calculate the sum of sqr sqrt
% 
% 
% 
% 
% 
% 
%   cutV=mean(Dat(BM)) -2*std(Dat(BM));
%     %noise estraction - we will need to genralize this somehow espesialy
% 
%     % selecting the relevant slices (this is not beutiful)
%     HM=ones(size(BM));
%     HM(Dat<cutV)=0;
%       % feeling the holes in the mask.
%     for i=1:size(HM,3)
%         HM(:,:,i)=imfill(HM(:,:,i),'holes');
%     end;
%    
% HM= HM & PDresamp>0;
% 
% 
%  G=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0','AllCoilsM0_phantomExample.nii.gz'));
% Gsos=sqrt(sum(G.data.^2,4)); % calculate the sum of sqr sqrt
% 
% 
% PDresamp=PDresamp(:,:,1:90);
% HM=HM(:,:,1:90);
% Gsos=Gsos(16:105,6:112,:);
% HM= HM & Gsos>0;
% for ii=1:size(G.data,4);
%     
%    Gi(:,:,:,ii)=G.data(16:105,6:112,ii);
% end
% 
% 
% 
% Gainfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain', 'new','SimValues','Gain.nii.gz');
% 
% dtiWriteNiftiWrapper(single(Gi),M0UnderSamp_Xform,Gainfile);
% 
% PDfile=fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain', 'new','SimValues','PD.nii.gz');
% dtiWriteNiftiWrapper(single(PDresamp),M0UnderSamp_Xform,Gainfile);
% 
% 
% 
