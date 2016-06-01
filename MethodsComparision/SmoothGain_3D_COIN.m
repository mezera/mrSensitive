function [PDfit1,Gfit]= SmoothGain_3D_COIN(M0file,maskfile,outputhPath,smoothKarnel)
% [PDfit1,Gfit]= SmoothGain_3D_COIN(M0file,maskfile)
% Use COIN code to estimate PD and coil GAIN form M0 image
% Noterdaeme, O., Anderson, M., Gleeson, F., & Brady, S. M. (2009). Intensity correction with a pair of spoiled gradient recalled echo images. Physics in Medicine and Biology, 54, 3473?3489. doi:10.1088/0031-9155/54/11/013
% we have modified the original pipeline :
%1. We have used a  mask to define the area that is approximately constant. The mask can be asegmented tissue like the white matter from the T1 map.
%2. Since COIN is 2D algorithm we have run it on the tree different plains  (axial, sagittal, coronal) and join all the 2D fits.
%
% Input: 
% M0file   - The path to the M0 image (3D) 
% Maskfile - The path to the mask image (3D). The mask is a region we
%           assume is about constant so that it?s smooth variance is the gain bias filed.
%outputhPath - location were to save the output as nifti (if empty the codewill run and won't save)
%smoothKarnel - a 3D smooth karnell of M0 data
% OutPut
% Gfit   - The fitted filed of M0
% PDfit1 - The fitted unbais PD map PDfit1=M0/Gfit.
%
% % Copyright Aviv Mezer , The Hebrew University 2014

%%
% load data 
PD=readFileNifti(M0file);xform=PD.qto_xyz; PD=PD.data;

if ~notDefined('smoothKarnel')
PD=smooth3(PD,'gaussian',smoothKarnel);
end
mask=readFileNifti(maskfile);mask=logical(mask.data);
BM=PD>0;
mask=mask & PD>0 ;
PD=PD./median(PD(mask));
Gi=zeros(size(PD));
Gi=PD;
Gi(~mask)=nan;
Gi(isinf(Gi))=nan;

%% Z slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));

Gz=zeros(size(PD));


%loop over slices
for  jj=1:size(PD,3)
    
    tmp=Gi(:,:,jj);
    tmpBM=BM(:,:,jj);
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if    (length(find(tmp>0))/length(find(tmpBM>0))>0.1  && length(wh)>80)
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth vertion of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gz(:,:,jj)=ZI;
        
        clear ZI
    end;
    
end


PDz=PD./Gz;
PDz(PDz>2)=0; PDz(PDz<0)=0;

%%  Y slice

% we will estimate  a smooth function in 2D space (Y -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));

Gy=zeros(size(PD));


%loop over slices
for  jj=1:size(PD,2)
    
    tmp=squeeze(Gi(:,jj,:));
    tmpBM=squeeze(BM(:,jj,:));
    %check that there is data in the slice
    wh=find(tmp>0);
    if   (length(find(tmp>0))/length(find(tmpBM>0))>0.1  && length(wh)>80);
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth vertion of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gy(:,jj,:)=ZI;
        
        clear ZI
    end;
    
end;


PDy=PD./Gy;

PDy(PDy>2)=0; PDy(PDy<0)=0;

%%  X slice

% we will estimate  a smooth function in 2D space (X -slice)
[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));

Gx=zeros(size(PD));


%loop over slices
for  jj=1:size(PD,1)
    
    tmp=squeeze(Gi(jj,:,:));
    tmpBM=squeeze(BM(jj,:,:));
    
    %check that there is data in the slice
    wh=find(tmp>0);
    if   (length(find(tmp>0))/length(find(tmpBM>0))>0.1  && length(wh)>80);
        %find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % estimate a smooth vertion of the data in the slice for original code see:
        % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        % put the result gain in the 3D gain image and fix orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gx(jj,:,:)=ZI;
        
        clear ZI
    end;
    
end;


PDx=PD./Gx;
PDx(PDx>2)=0; PDx(PDx<0)=0;

%% join the three dimention PD
% We will average the position that PD have estimations.

PDfit=zeros(size(PD));
PDvox=zeros(size(PD));

whx=find(PDx>0);
PDfit(whx)=PDx(whx);
PDvox(whx)=PDvox(whx)+1;


why=find(PDy>0);
PDfit(why)=PDfit(why)+PDy(why);
PDvox(why)=PDvox(why)+1;



whz=find(PDz>0);
PDfit(whz)=PDfit(whz)+PDz(whz);
PDvox(whz)=PDvox(whz)+1;

PDfit=PDfit./PDvox;


G=PD./PDfit;

G(isinf(G))=nan;
G(G<0)=nan;

Gfit=zeros(size(PD));
mask=PDfit>0 & BM;

%% estimate G
% After running on Z Y X and average the G fit (above).
% We will run on G one more time to assure smooth and complete G fit.
% We will do it again on each of the three plains

%% Gz
% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));


    %loop over slices
    for  jj=1:size(PD,3)
        
        tmp=G(:,:,jj);
        tmpBM=mask(:,:,jj);
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>100);
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
            % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            % put the result gain in the 3D gain image and fix orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
  
%% Gx
% we will estimate  a smooth function in 2D space (X -slice)
[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));



    for  jj=1:size(PD,1)
        
        tmp=squeeze(Gfit(jj,:,:));
        
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>80);
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
            % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            % put the result gain in the 3D gain image and fix orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
    


%% Gy

[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));

%loop over slices
    for  jj=1:size(PD,2)
        
        tmp=squeeze(Gfit(:,jj,:));
        
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>80);
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
            % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            % put the result gain in the 3D gain image and fix orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(:,jj,:)=ZI;
            
            clear ZI
        end;
        
    end;
   
Gfit(isinf(Gfit))=0;
Gfit(Gfit<0)=0;
Gfit(isnan(Gfit))=0;



%% having a final Gain function  we will calculate PD

PDfit1=PD./Gfit;

PDfit1(isinf(PDfit1))=0;
PDfit1(PDfit1<0)=0;
PDfit1(isnan(PDfit1))=0;
PDfit1(PDfit1>4)=4;
%%
if ~notDefined('outputhPath')
    PDfile=fullfile(outputhPath,'PDfit.nii.gz');
        dtiWriteNiftiWrapper(single(PDfit1),xform,PDfile);
      Gfile=fullfile(outputhPath,'Gfit.nii.gz');
        dtiWriteNiftiWrapper(single(Gfit),xform,Gfile);
end

    