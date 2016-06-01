function [opt, G, PD]=mrQ_smoothGain_step4a(opt,PD)

%multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
%T1

%Brain mask
BM=readFileNifti(opt.BMfile);
xform=BM.qto_xyz;
BM=BM.data;

%%
%calculate coil gains from the PD_fit and M0. M0 is 4D image the 4th
%dimention is the image of each coil
% M0=Gain X PD  --> Gain= M0 / PD;
% PD fit is not full image. we have values only where the boxes fit where done
% and got a "good sulotion". hopfully we have most of the image.
% It is given that the gain function vart smoothly in space.
% Therefore to fill the area where PD data is missing we will model the
% Gain as a smooth function in space.


% This model will also smooth area where the gain estimation is noise.
mask=BM & PD>0 ;


%% Z slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));

G=zeros(size(M0));


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
    
    
    %the estimate gain for the coil
    Gi=zeros(size(PD));
    Gi=M0(:,:,:,ii)./PD;
    Gi(~mask)=nan;
    Gi(isinf(Gi))=nan;
    
    
    %loop over slices
    for  jj=1:size(PD,3)
        
        tmp=Gi(:,:,jj);
        tmpBM=BM(:,:,jj);
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if    (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800)
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
            Gi(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end


PDz=M0./G;
PDz=median(PDz,4);

Gz=G;

%%  Y slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));

G=zeros(size(M0));


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
    
    
    %the estimate gain for the coil
    Gi=zeros(size(PD));
    Gi=M0(:,:,:,ii)./PD;
    Gi(~mask)=nan;
    Gi(isinf(Gi))=nan;
    
    
    %loop over slices
    for  jj=1:size(PD,2)
        
        tmp=squeeze(Gi(:,jj,:));
        tmpBM=squeeze(BM(:,jj,:));
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800);
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
            Gi(:,jj,:)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end


PDy=M0./G;
PDy=median(PDy,4);
Gy=G;


%%  X slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));

G=zeros(size(M0));


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
    
    
    %the estimate gain for the coil
    Gi=zeros(size(PD));
    Gi=M0(:,:,:,ii)./PD;
    Gi(~mask)=nan;
    Gi(isinf(Gi))=nan;
    
    
    %loop over slices
    for  jj=1:size(PD,1)
        
        tmp=squeeze(Gi(jj,:,:));
        tmpBM=squeeze(BM(jj,:,:));
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800);
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
            Gi(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end


PDx=M0./G;
PDx=median(PDx,4);
Gx=G;

 %% join the three dimention PD
% 
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









%% join the 3 dimention of Gain estimation 
% 
% Gfit=zeros(size(M0));
% Gvox=zeros(size(M0));
% 
% whx=find(Gx>0);
% Gfit(whx)=Gx(whx);
% Gvox(whx)=Gvox(whx)+1;
% 
% why=find(Gy>0);
% Gfit(why)=Gfit(why)+Gy(why);
% Gvox(why)=Gvox(why)+1;
% 
% whz=find(Gz>0);
% Gfit(whz)=Gfit(whz)+Gz(whz);
% Gvox(whz)=Gvox(whz)+1;
% 
% Gfit=Gfit./Gvox;
% 
% PDfit1=M0./Gfit;
% PDfit1=median(PDfit1,4);
%%

%% estimate G
%% Gz

mask=PDfit>0 & BM;
% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));

G=zeros(size(M0));


%l run over coils and fit the Gain for each

for ii=1:size(M0,4)
    
    
    %the estimate gain for the coil
    Gi=zeros(size(PD));
    Gi=M0(:,:,:,ii)./PDfit;
    Gi(~mask)=nan;
    Gi(isinf(Gi))=nan;
    
    
    %loop over slices
    for  jj=1:size(PD,3)
        
        tmp=Gi(:,:,jj);
        tmpBM=mask(:,:,jj);
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>1000);
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
            Gi(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end



%% Gx
% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
    clear Gi
    Gi=G(:,:,:,ii);
    %loop over slices
    for  jj=1:size(PD,1)
        
        tmp=squeeze(Gi(jj,:,:));
        
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>800);
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
            Gi(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end

%% Gy

[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
     clear Gi
    Gi=G(:,:,:,ii); 
    %loop over slices
    for  jj=1:size(PD,2)
        
        tmp=squeeze(G(:,jj,:));
        
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>800);
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
            Gi(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end


%% and back for PD
 PDfit1=M0./G;
PDfit1=median(PDfit1,4);
%% save PD and Coil gain & Upsample if needed to the original resulotion

% save the gain in the resulotion we calculate it
G_filename=fullfile(opt.outDir, 'Gains.nii.gz');
dtiWriteNiftiWrapper(single(G),xform,G_filename);

PD_filename=fullfile(opt.outDir, 'PD.nii.gz');


if ~isfield(opt,'Resamp')
    opt.Resamp=0;
end
if ~opt.Resamp==1; %no resample we will save PD and coil gain
    dtiWriteNiftiWrapper(single(PDfit1),xform,PD_filename);
    
else
    %we will up sample the G to the original M0 size and then divied M0/G
    %to get PD
    sz=size(G);
    bb = mrAnatXformCoords(xform, [1 1 1; sz(1:3)]);
    clear M0
    % get the original M0
    M0=readFileNifti(opt.M0file_Org);outMm=M0.pixdim(1:3);     M0=M0.data; MSZ=size(M0);
    % loop over coils and up sample G
    for ii=1:sz(4)
        
        [Gi, M0UnderSamp_Xform] = mrAnatResliceSpm(double(G(:,:,:,ii)), inv(xform), bb, outMm, 1, 0);
        GSZ=size(Gi);
        % get the size difference is exsist
        Sizes=[min(GSZ(1),MSZ(1)) min(GSZ(2),MSZ(2)) min(GSZ(3),MSZ(3)) ]; % we calculate the sizes of the two image M0 and G becouse resample may end up with adding or losing a line of zeros.
        
        % calculate the coil PD
        M0(1:Sizes(1),1:Sizes(2),1:Sizes(3),ii)=M0(1:Sizes(1),1:Sizes(2),1:Sizes(3),ii)./Gi(1:Sizes(1),1:Sizes(2),1:Sizes(3));
        
    end
    PD=median(M0,4); %is wight sum by SNR better
    
    dtiWriteNiftiWrapper(single(PD),M0UnderSamp_Xform,PD_filename);
    
end

opt.PDfile=PD_filename;
opt.Gainfile=G_filename;
save(opt.logname,'opt')
