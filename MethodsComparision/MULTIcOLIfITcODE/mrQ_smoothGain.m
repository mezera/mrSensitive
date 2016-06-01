function [Gfit]=mrQ_smoothGain(Gi,mask)
%%
%calculate smooth coil gains from the gain estimation. 
% Gain fit is not full image (i.e. only mask). 
% It is given that the gain function vart smoothly in space.
% Therefore to fill the area where data is missing we will model the
% Gain as a smooth function in space.


% This model will also smooth area where the gain estimation is noise.

Gi(~mask)=nan;
Gi(isinf(Gi))=nan;

%% Z slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(Gi,1),1:size(Gi,2));

Gz=zeros(size(Gi));


%l run over coils and fit the Gain for each


%the estimate gain for the coil


%loop over slices
for  jj=1:size(Gi,3)
    
    tmp=Gi(:,:,jj);
    tmpBM=mask(:,:,jj);
    
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
        Gz(:,:,jj)=ZI;
        
        clear ZI
    end;
    
    
    
end


%%  Y slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(Gi,1),1:size(Gi,3));

Gy=zeros(size(Gi));

%the estimate gain for the coil



%loop over slices
for  jj=1:size(Gi,2)
    
    tmp=squeeze(Gi(:,jj,:));
    tmpBM=squeeze(mask(:,jj,:));
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
        Gy(:,jj,:)=ZI;
        
        clear ZI
    end;
    
end;





%%  X slice

% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(Gi,2),1:size(Gi,3));

Gx=zeros(size(Gi));


%loop over slices
for  jj=1:size(Gi,1)
    
    tmp=squeeze(Gi(jj,:,:));
    tmpBM=squeeze(mask(jj,:,:));
    
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
        Gx(jj,:,:)=ZI;
        
        clear ZI
    end;
    
end;



%% join the three dimention PD
%
Hval=prctile(Gi(mask),99);
Lval=prctile(Gi(mask),1);

Gfit=zeros(size(Gi));
Gvox=zeros(size(Gi));

whx=find(Gx>Lval & Gx<Hval);
Gfit(whx)=Gx(whx);
Gvox(whx)=Gvox(whx)+1;

why=find(Gy>Lval & Gy<Hval);
Gfit(why)=Gfit(why)+Gy(why);
Gvox(why)=Gvox(why)+1;


whz=find(Gz>Lval & Gz<Hval);
Gfit(whz)=Gfit(whz)+Gz(whz);
Gvox(whz)=Gvox(whz)+1;

Gfit=Gfit./Gvox;


Gfit(isinf(Gfit))=nan;
Gfit(Gfit<0)=nan;

mask=Gfit>0 & mask;

%% estimate G
%% Gz
% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(Gi,1),1:size(Gi,2));


   
    
    %loop over slices
    for  jj=1:size(Gi,3)
        
        tmp=Gfit(:,:,jj);
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
            Gfit(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
   


%% Gx
% we will estimate  a smooth function in 2D space (Z -slice)
[XI, YI]=meshgrid(1:size(Gi,2),1:size(Gi,3));



    for  jj=1:size(Gi,1)
        
        tmp=squeeze(Gfit(jj,:,:));
        
        
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
            Gfit(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;
    
    
    

%% Gy

[XI, YI]=meshgrid(1:size(Gi,1),1:size(Gi,3));


%l run over coils and fit the Gain for each
%loop over slices
    for  jj=1:size(Gi,2)
        
        tmp=squeeze(Gfit(:,jj,:));
        
        
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
            Gfit(:,jj,:)=ZI;
            
            clear ZI
        end;
        
    end;
   
Gfit(isinf(Gfit))=0;
Gfit(Gfit<0)=0;
Gfit(isnan(Gfit))=0;


