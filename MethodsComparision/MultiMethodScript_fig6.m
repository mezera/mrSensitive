%  M0 multi channels local fitting methods for figure 6 in Mezer et. al.
%  HBM 2016

% AM/BW Vistaosft Team, 2013

%% Calculate on PD on  a largre phantm voulume
% running the solotion  for PD from raw image will take more then 48
% hours on single CPU and few hours with SGE paralell computing. To stat
% form Raw please see : https://github.com/mezera/mrQ

disp('code for simulations in figure 6')


%%2. fit with T1 regularization (change starting point) change number of
%%coils 4,8 32

%% 3. fit with corralation loop over 1.starting point), 2.number of coils (4,8,32) 3. number of tissue (1:4) 5. fitting method (T1 reg, Corr ,ridge)

%% make a mrW stracture:

%% to solve a full voulume we will use mrQ software
% soloving this problem can be very slow (~24 hours on a single CPU) with
% out SGE parallel computing.

% We show here the script form the article that use SOS coil summation and for PD free of true smooth PD in space.
% to plot to Median summationchange to Median when ever it say SOS and 
% to change to smooth PD in space change 'new' whenever we use'old'

% make a mrQ directory and parameters
M0file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old' ,'InPut','M0_noise.nii.gz'));
%1. make M0 with G PD and noise have T1
T1file=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old','InPut','T1.nii.gz'));
BMfile=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain' ,'old', 'InPut','mask.nii.gz'));

dirpath =*** Path_for_outPuts ***



%% set the inputs form the simulated PhantomBrain

%% Intiate the multi box fit parameters

for Reg=1:3 % differnt methods of regularization 1-T1  2-corr 3-Tikonov. See in mrQ_CoilPD_gridFit_Multi.m
    for Ncoils=1:1
        
        if Ncoils==1; Coilsinfo.maxCoil=4;Coilsinfo.minCoil=4; Coilsinfo.useCoil=[1:16]; end % use data of 4 out of 16 coils
        if Ncoils==2; Coilsinfo.maxCoil=8;Coilsinfo.minCoil=8; Coilsinfo.useCoil=[1:16]; end% use data of 8 out of 16 coils
        if Ncoils==3; Coilsinfo.maxCoil=32;Coilsinfo.minCoil=32; Coilsinfo.useCoil=[1:32]; end% use data of 32 out of 32 coils
        
        for Init=3:3 % different intiations are optional
            
            run_mrQtest_PDfit(dirpath,M0file,BMfile,T1file,Coilsinfo,Reg,Init)
        end
    end
end


%% build the PD map
dirpath ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/PDFItMethodTesting_multiCondisions'

%dirpath ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting'
for Reg=3:3
    if Reg==1 tresh=0.01; CutResid=0.1; else, tresh=0.05; CutResid=0.95; end
    
    for Ncoils=1:1
        
        if Ncoils==1; Coilsinfo.maxCoil=4;Coilsinfo.minCoil=4; Coilsinfo.useCoil=[1:16]; end
        if Ncoils==2; Coilsinfo.maxCoil=8;Coilsinfo.minCoil=8; Coilsinfo.useCoil=[1:16]; end
        if Ncoils==3; Coilsinfo.maxCoil=32;Coilsinfo.minCoil=32; Coilsinfo.useCoil=[1:32]; end
        
        for Init=3:3
            
            mrQ_file=  fullfile(dirpath,['mrQ_R_' num2str(Reg) '_C_' num2str(Coilsinfo.maxCoil) '_I_' num2str(Init)],'mrQ_params.mat')
            load (mrQ_file)
            mrQ.opt=mrQ_buildPD_ver2(mrQ.opt_logname,0,[],CutResid,[],tresh);
            mrQ.SPGR_PDBuild_done=1;
            save(mrQ.name,'mrQ');
        end
    end
end



%% T1_PD build single coil
%% T1 PD single coil
Coilsinfo.maxCoil=1;Coilsinfo.minCoil=1; Coilsinfo.useCoil=1;
Reg=5;Init=5; %i
dirpath1 ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/PDFItMethod_OnlyT1'
%dirpath ='/biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting'
M0file=fullfile(dirpath1,'M0');
run_mrQtest_PDfit(dirpath,M0file,BMfile,T1file,Coilsinfo,Reg,Init)
  

load  /biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/PDFItMethodTesting_multiCondisions/mrQ_R_5_C_1_I_5/mrQ_params.mat

            mrQ.opt=mrQ_buildPD_ver2(mrQ.opt_logname,0,[],[],[],0.01);
 %new         
 load  /biac4/wandell/biac2/wandell2/data/WMDevo/code/Sensitive/PDFitMethodTesting/mrQ_R_5_C_1_I_6/mrQ_params.mat
 
  mrQ.opt=mrQ_buildPD_ver2(mrQ.opt_logname,0,[],[],[],0.01);
 

