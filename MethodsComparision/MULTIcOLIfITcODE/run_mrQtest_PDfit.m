function mrQ=run_mrQtest_PDfit(dir,M0file,BMfile,T1file,Coilsinfo,Reg,Init,SunGrid,proclus)

if notDefined('dir')
dir=(fullfile(mrSensitiveRootPath,'ExampleData','PhantomBrain') )
end


if notDefined('Coilsinfo')
Coilsinfo.maxCoil=4;
    Coilsinfo.minCoil=4;
    Coilsinfo.useCoil=[1:16];
end

if notDefined('Reg')
Reg=1; % Reg=1 T1reg, Reg=0 corralation ; Reg=2 ridge
end


if notDefined('Init')
        Init=1;
end

if notDefined('StartingMethod')
StartingMethod='sos';
end

%% make a mrQ stracture
 outDir = fullfile(dir,['mrQ_R_' num2str(Reg) '_C_' num2str(Coilsinfo.maxCoil) '_I_' num2str(Init)]);

            if ~exist(outDir,'dir'); mkdir(outDir); end
            mrQ = mrQ_Create(dir,[],outDir);
            
% Let's skip all an relevant mrQ processing  until the PD and coil
% sensativity mapping
            mrQ.Arange_Date=date;
            mrQ.SEIR_done=1;
            mrQ.SPGR_init_done=1;
            mrQ.SPGR_coilWeight_done=1;
            mrQ.SPGR_T1fit_done=1;
            mrQ.segmentaion=1;
            mrQ.calM0_done=1;
            mrQ.calM0_done=1;
            mrQ.brakeAfterPD=1;

% let save the M0 image
mrQ.M0combineFile=M0file;


%% if using SGE 
if notDefined('SunGrid')
mrQ.SunGrid=1;
else
    mrQ.SunGrid=SunGrid;
end
%if not  not evaluate:
%mrQ.SunGrid=0;

%% if using the STANFORD proclus computer 

if notDefined('proclus')
mrQ.proclus=1;
else
    mrQ.proclus=proclus;
end
%if not evaluate:
%mrQ.proclus=0;

% where we  read and write to
mrQ.spgr_initDir=outDir;

save(mrQ.name,'mrQ');


%%






[mrQ.opt_logname]=mrQ_PD_multicoil_RgXv_GridCall_Fittesting(mrQ.spgr_initDir,mrQ.sub,mrQ.PolyDeg,M0file,T1file,BMfile,[],[],[],Coilsinfo,Reg,Init,[],mrQ);
%
save(mrQ.name,'mrQ');
%% solve for each box (via SGE or not)
%mrQ_fitM0boxesCall(mrQ.opt_logname,mrQ.SunGrid,mrQ.proclus);
%mrQ_fitM0boxesCall_ver3(mrQ.opt_logname,mrQ.SunGrid,mrQ.proclus);
%for only local T1 no multi coils
if Reg==5
mrQ_fitM0boxesCall_T1PD(mrQ.opt_logname)
else
% for multi coil with different regularizations
mrQ_fitM0boxesCall_MultiTest(mrQ.opt_logname,mrQ.SunGrid,mrQ.proclus);
end
return
%% build the from the results
%mrQ_run(mrQ.name);
%mrQ_run_Ver1(mrQ.name);