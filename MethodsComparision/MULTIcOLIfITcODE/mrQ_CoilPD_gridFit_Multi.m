function mrQ_CoilPD_gridFit_Multi(opt,jumpindex,jobindex)
%
% mrQ_CoilPD_gridFit_ver3(opt,jumpindex,jobindex)
%  this function call by the sun grid it load the relavant data and fit the
%  PD and coils bias of M0 mage rigion (voulume).
% the imaging voulume region also call here "box". The box is a location (few voxel 100's to
% 1000's).  the idea is to use the information in each coils. we tery to find
% to PD that is similar to all the coils images and fit the coil bias that
% is diferent for each coil.
%  The fit is done in few steps.
%  1. we get the box data.
% 2. we fit subset of coils informatio nwith regularizartion (of T1) and crossvalidation to it's quallety
% 4. the best cross validation regularization is used to fit all the data.
%4. we save the resuts
%
% INPUTS:
%       opt - this is optmization structure that was passed from
%       mrQ_fitPD_multicoil and it have all the needed information
%       jumpindex - how many boxes this grid call we fit (book keeping)
%       jobindex  - the number of box it will start whe ncalling the grid
%       (book keeping)
%
% OUTPUTS:
%  save an output file with fitted parameters in a tmp directorry
%   this will be used lster by mrQfitPD_multiCoils_M0 to make the PD map

% SEE ALSO:
% mrQ_PD_multicoil_RgXv_GridCall
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization



%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%cheack that this box have brain data
if ed>length(opt.wh), ed=length(opt.wh);end;

nIteration=ed-st+1;
%intilazie parameters and saved outputs

% Get the M0 and T1 information

%multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
%T1
T1=readFileNifti(opt.T1file);
T1=T1.data;


%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

%seg mask
seg=readFileNifti(opt.segfile);
seg=seg.data;


smoothkernel=opt.smoothkernel;


% thepoly basis to fit the coil gains
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);
maxCoil=opt.maxCoil;
minCoil=opt.minCoil;
useCoil=opt.useCoil;
nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

% initite the saved parameters
fb=zeros(nIteration,1,3);
gEst=zeros(nPolyCoef,maxCoil,nIteration);
resnorm=zeros(nIteration,1);
exitflag=zeros(nIteration,1);
skip=zeros(nIteration,1);
X_valdationErrSN=zeros(nIteration,1);
if isfield(opt, 'lambda')
    X_valdationErr=zeros(2,length(opt.lambda),nIteration);
else
    X_valdationErr=[];
end


BestReg=zeros(nIteration,2);
Clists=zeros(maxCoil, nIteration);
ResidErr=zeros(nIteration,1);
Clists2=zeros(maxCoil, nIteration);
BiasTest=0;
Iter=0;

%%  II. go over the box the boxs


for ii= st:ed,
    %run over the box you like to fit
    clear M01  t1  BM1  SZ M0_v R1basis PDinit Segmask g0 G0 mask1 X_valdationErrF   X_gEstF  best2 Segmask
    Iter= Iter+1;
    tic
    %find the x,y,z location of the box (this is not x,y,z location in image space but
    %grid of boxes we made by mashgrid in  mrQ_PD_multicoil_RgXv_GridCall.m
    [fb(Iter,1,1), fb(Iter,1,2), fb(Iter,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, t1, BM1, SZ, skip(Iter), Segmask]= mrQ_GetM0_boxData(opt,T1,M0,BM,fb(Iter,1,:),smoothkernel,seg);
    M0_v = reshape(M01, prod(SZ(1:3)), SZ(4));
    R1=1./(t1(:));
    
    if  skip(Iter)==1
        disp(['skipping box ' num2str(ii) ' bad data'])
        
    else
        
        
        
        
        %% Select the coils to fit
        %
        if length(useCoil)>SZ(4); useCoil=useCoil(1:SZ(4));end
        
        Clist=mrQ_select_coils(maxCoil,max(useCoil),M0_v(BM1,:));
        nCoils=length(Clist);
        Clists(1:nCoils,Iter)=Clist;
        
        % find an alternative coil list
        if nCoils<=SZ(4)/2
            BiasTest=true;
            
            Mtmp=M0_v(BM1,:);
            Mtmp(:,Clist)=1;
            Clist2=mrQ_select_coils(maxCoil,max(useCoil),Mtmp);
            Clists2(1:nCoils,Iter)=Clist2;
            clear Mtmp;
        else
            BiasTest=false;
        end
        
        
        % find an alternative coil list
        if nCoils<SZ(4)/3
            BiasReTest=true;
            
            Mtmp=M0_v(BM1,:);
            Mtmp(:,Clist)=1;
            Mtmp(:,Clist2)=1;
            Clist3=mrQ_select_coils(maxCoil,max(useCoil),Mtmp);
            Clists3(1:nCoils,Iter)=Clist3;
            clear Mtmp;
        else
            BiasReTest=false;
        end
        
        %%
        
        
        % The nonlinsqr fit with ridge fits
        [ PDinit,g0] =Get_PDFit_InIt(opt.Init,M0_v(:,Clist),pBasis,R1,Segmask)      ;
        
        options = optimset('Display','off',...  %'iter'final
            'MaxFunEvals',Inf,...
            'MaxIter',200,...
            'TolFun', 1e-6,...
            'TolX', 1e-6,...
            'Algorithm','levenberg-marquardt');
        
        
        
        
        
        %% Select the method to Fit
        
        if opt.Reg==1
            %%  X-Validation Fit of coil Gain with T1 regularization
            %make a R1 regularazation matrix
            R1basis(1:nVoxels,1) = 1;  R1basis(:,2) =R1;
            R1basis=double(R1basis);
            
            [X_valdationErr(:,:,Iter), BestReg(Iter,1), BestReg(Iter,2),X_valdationErrSN(Iter), PDfit, gEst(:,:,Iter), resnorm(Iter),exitflag(Iter), G] = ...
                mrQ_FitPDT1regXCall(opt,M0_v,BM1,pBasis,Clist,R1basis,g0,Segmask);
            
            
            if BiasTest
                % get ht ePD with a different set of coils
                [PDfit2 ]=pdCoilSearch_T1reg( opt.lambda(BestReg(Iter,2)),M0_v(BM1,Clist2),pBasis(BM1,:),R1basis(BM1,:),[],Segmask(BM1));
                
                %%
                tmp=zeros(SZ(1:3));
                tmp(BM1)=PDfit2(:);
                PDfit2=tmp;
                tmp(BM1)=PDfit(:);
                PDfit=tmp;
                
                PDfit2  = reshape(PDfit2,SZ(1:3));
                PDfit  = reshape(PDfit,SZ(1:3));
                
                ErrMap=PDfit./mean(PDfit(BM1))-PDfit2./mean(PDfit2(BM1));
                [~,~,ResidErr(Iter)] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);
            end
            if BiasReTest & ResidErr(Iter)>0.1
                % in case the coil does not agree (big bias) and we can find a third set of coils we will try to find sulotion that agrea.
                [ PDinit,g0] =Get_PDFit_InIt(opt.Init,M0_v(:,Clist2),pBasis,R1,Segmask)      ;
                [X_vald2, BestReg2(1), BestReg2(2),X_valdationErrSN2, PDfit2, gEst2, resnorm2,exitflag2, G] = ...
                    mrQ_FitPDT1regXCall(opt,M0_v,BM1,pBasis,Clist2,R1basis,g0,Segmask);
                
                [ PDinit,g0] =Get_PDFit_InIt(opt.Init,M0_v(:,Clist3),pBasis,R1,Segmask)      ;
                [X_vald3, BestReg3(1), BestReg3(2),X_valdationErrSN3, PDfit3, gEst3, resnorm3,exitflag3, G] = ...
                    mrQ_FitPDT1regXCall(opt,M0_v,BM1,pBasis,Clist3,R1basis,g0,Segmask);
                tmp=zeros(SZ(1:3));
                tmp(BM1)=PDfit2(:);
                PDfit2=tmp;
                PDfit2  = reshape(PDfit2,SZ(1:3));
                
                
                 tmp=zeros(SZ(1:3));
                tmp(BM1)=PDfit3(:);
                PDfit3=tmp;
                PDfit3  = reshape(PDfit3,SZ(1:3));
                
                % now we can check if a new set of coils data still end up with a bias
                
                ErrMap31=PDfit./mean(PDfit(:))-PDfit3./mean(PDfit3(:));
                [~,~,ResidErr31] = fit3dpolynomialmodel(ErrMap31,logical(ErrMap31),1);
                
                ErrMap23=PDfit3./mean(PDfit3(:))-PDfit2./mean(PDfit2(:));
                
                [~,~,ResidErr23] = fit3dpolynomialmodel(ErrMap23,logical(ErrMap23),1);
                
                if ResidErr31<0.1 % if 1 and 3 agree then 2 was probably worng
                    %              we need to save the new coils sets
                    ResidErr(Iter)=ResidErr31;
                    Clist2=Clist3;
                    Clists2(1:nCoils,Iter)=Clist3;
                    
                elseif ResidErr23<0.1 & ResidErr31>0.1    %if 2 3 agrea then 1 was probably wrong
                    %             Then we need to save the new coils sets and
                    %             parameters (replacing set no 1
                    
                    X_valdationErr(:,:,Iter)=X_vald2; BestReg(Iter,:)=[BestReg2(1), BestReg2(2)];
                    gEst(:,:,Iter)=gEst2; X_valdationErrSN(Iter)=X_valdationErrSN2;
                    resnorm(Iter)=resnorm2; exitflag(Iter)=exitflag2;
                    
                    Clists(1:nCoils,Iter)=Clist2;
                    Clists2(1:nCoils,Iter)=Clist3;
                    
                end
            end
        end
        
        
        
        
        
        
        %% Use Coil correlation regularization
        if opt.Reg==2;
            coefdat = tril(corrcoef(M0_v(BM1,Clist)),-1);
            
            RegWeight  = 1000;
            TissueMask = logical(M0_v(:,1));
            XvalidationMask = logical(M0_v(BM1,Clist));
            nUseCoils=length(Clist)
            [gEst(:,:,Iter), resnorm(Iter),~,exitflag(Iter)] = ...
                lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,double(M0_v(:,Clist)),pBasis,nUseCoils,RegWeight,BM1,double(coefdat),XvalidationMask),g0,[],[],options);
            
            [PDfit, Gn] = pdEstimate(M0_v(:,Clist),pBasis, gEst(:,:,Iter));
            
            
            if BiasTest
                % get ht ePD with a different set of coils
                coefdat = tril(corrcoef(M0_v(BM1,Clist2)),-1);
                
                [g, resnorm1,dd1,exitflag1] = ...
                    lsqnonlin(@(par) errFitNestBiLinearCorrReg(par,double(M0_v(:,Clist2)),pBasis,nUseCoils,RegWeight,BM1,double(coefdat),XvalidationMask),gEst(:,:,Iter),[],[],options);
                
                [PDfit2, Gn] = pdEstimate(M0_v(:,Clist2),pBasis, g);
                
                %%
                tmp=zeros(SZ(1:3));
                tmp(BM1)=PDfit2(BM1);
                PDfit2=tmp;
                tmp(BM1)=PDfit(BM1);
                PDfit=tmp;
                
                PDfit2  = reshape(PDfit2,SZ(1:3));
                PDfit  = reshape(PDfit,SZ(1:3));
                
                ErrMap=PDfit./mean(PDfit(:))-PDfit2./mean(PDfit2(:));
                [~,~,ResidErr(Iter)] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);
                
                
                
            end
            
        end
        
        
        if opt.Reg==3;
            maxLoops = 1000;
            sCriterion = 1e-4;  % Stopping criterion
            
            % Set the fiiting loop parmeters
            [X_valdationErrF,  ]=pdX_valdationLoop_RidgeReg_ver1(opt.lambda,3,M0_v(BM1,Clist), pBasis(BM1,:),PDinit(BM1),maxLoops,sCriterion);
            X_valdationErr(:,:,Iter)=X_valdationErrF;
            %mrvNewGraphWin;loglog(opt.lambda+1e-1,X_valdationErr(2,:,Iter),'*-'); xlabel('lambda');ylabel('X-V error');
            best1 = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:))); % sum of abs err
            best2 = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)));% RMSE
            
            BestReg(Iter,:)=[best1(1) best2(1)];
            
            
            
            
            RidgeBLFit = pdBiLinearRidgeFit(M0_v(BM1,Clist), pBasis(BM1,:),opt.lambda(best2),maxLoops,sCriterion,PDinit(BM1));
            gEst(:,:,Iter)=RidgeBLFit.g;
            
            if BiasTest
                RidgeBLFit2 = pdBiLinearRidgeFit(M0_v(BM1,Clist2), pBasis(BM1,:),opt.lambda(best2),maxLoops,sCriterion,PDinit(BM1));
                tmp=zeros(SZ(1:3));
                tmp(BM1)=RidgeBLFit2.PD;
                PDfit2=tmp;
                tmp(BM1)=RidgeBLFit.PD;
                PDfit=tmp;
                
                PDfit2  = reshape(PDfit2,SZ(1:3));
                PDfit  = reshape(PDfit,SZ(1:3));
                
                ErrMap=PDfit./mean(PDfit(:))-PDfit2./mean(PDfit2(:));
                [~,~,ResidErr(Iter)] = fit3dpolynomialmodel(ErrMap,logical(ErrMap),1);
                
                
            end
        end
        
        
        toc
    end
end

name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'gEst','resnorm','exitflag','st','ed','skip','fb' ,'X_valdationErr','X_valdationErrSN','BestReg','Clists','Clists2','ResidErr' ,'BiasTest')


