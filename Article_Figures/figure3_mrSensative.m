%% figure 3
clear
disp('Figure 3')


%  hougenus phantom M0 data fitted with polynomials.
% Code associated with Mezer, et. al. 2016, HBM
%
% AM/BW Mezer Lab & Vistaosft Team, 2013


%%  Make sure mrSensitive is on the path
addpath(genpath(fullfile(mrSensitiveRootPath)));


%% Set parameters for the coils from the phantom data

nCoils   = 32;     % A whole bunch of coils
nDims    = 3;      % XYZ
noiseFloor = 500;  % This is the smallest level we consider
printImages  = false;   % No printing now
smoothkernel = [];      % Fit to the unsmoothed M0 data. NOt now
BasisFlag    = 'qr';    % Which matrix decomposition for fitting: Ortonormal

mrvNewGraphWin('Figure 3','wide');

%% Get the polynomial error
% loop over 5 (sampleLocation)  location along the phantom for edge to center and the other edge.
%loop over box voulume and polynomial order.
%in each case we fit the polynomial to each coil phantom M0 data
% we calculate the precent error btween the data and the  fit,


for sampleLocation = 6:10;% Which box
    %Loop over locations (the location and the phantom data are pre defined in phantomGetData.m)
    clear phantomP
    
    for nSamples=2:10
        %  Loop over box voulume The box is [-nSamples:nSamples -nSamples:nSamples -nSamples:nSamples] center around the location (sampleLocation)
        
        for pOrder=1:3
            % Loop over polynomial order. those are 3D ortonormal
            % polynomial define in polyCreateMatrix.m
            
            phantomP(nSamples,pOrder) = ...
                pdPolyPhantomOrder(nSamples, nCoils, nDims, pOrder, ...
                noiseFloor, sampleLocation, printImages, smoothkernel, BasisFlag);
            % pdPolyPhantomOrder.m get the data and polynomial fit the phantom data and svve the error. The result are save in a stacture phantomP
            
        end
    end


%% Extract key values from the structure containing the polynomial fits

% Number of samples is the box size
% pOrder is the polynomial order
PE = zeros(10,30);
Volume = zeros(10,3);

for nSamples=2:10 % The box is -nSamples:nSamples
    for pOrder=1:3 %  polynomial order
        
        % The box we use is a volume with -nSamples:nSamples on a side.
        % The resolution of the phantom scan voxel is 2mm. 
        % So multply by 2 and ^3 for volume.
        % The volume in mm3 is 
        Volume(nSamples,pOrder) = ((nSamples*2+1)*2)^3; 
% the mean percentError for a sample giving the order and voulue
        PE(nSamples,pOrder)     =  phantomP(nSamples,pOrder).percentError;

    end
end

%% plot the location precent error as function of voulume and polynomial order 
% Show the fit
N=sampleLocation-5
subplot(1,5,N)
hold on

plot(Volume(2:10,1),PE(2:10,1) ,'-k*', 'MarkerSize',10);

plot(Volume(2:10,2),PE(2:10,2) ,'-ko', 'MarkerSize',10);

plot(Volume(2:10,3),PE(2:10,3) ,'-ks', 'MarkerSize',10);
% loglog(Volume(2:10,3),ones(9,1),'--k')
if N==3
xlabel(' Volume mm^3','FontSize',20);
elseif N==1
    ylabel('percent error','FontSize',20);
end
%set(gca,'xscale','log','yscale','log','FontSize',24)
set(gca,'xscale','log','yscale','linear','FontSize',18)
%axis image; axis square
%title([' location no :'  num2str(sampleLocation)])
grid on
ylim([0 20])
set(gca,'YTick',[0 4 8 12 16 ])
set(gca,'YTickLabel',{'0' '4' '8' '12' '16'})

end

legend('1st order' , '2nd order', '3rd order','Location','NorthEast')


%% Ilustrate the location on sum of squre M0 image

M0=readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0','AllCoilsM0_phantomExample.nii.gz'));
M0=sqrt(sum(M0.data.^2,4)); % calculate the sum of sqr sqrt

% get a mask for area that are not outside the phantom or on the edge of
% it.
Mask =readFileNifti(fullfile(mrSensitiveRootPath,'ExampleData','AgarPhantom' ,'M0','PhantomMask.nii.gz'));
Mask=logical(Mask.data);


M0(~Mask)=0;

% mark the location as a black spots

 %the pre define location along the Y axes up to dwon see phantomGetData.m

   
M0(30, 60,45)=0; %6
M0(45, 60,45)=0; % 7
M0(60, 60,45)=0;% 8
M0(75, 60,45)=0;% 9
M0(90, 60,45)=0;% 10

%plot
showMontage(rot90(M0(:,:,45)));
axis square;axis off; colormap hot;colorbar off

%% gloobal fit
Mask= Mask & M0>0;
V=zeros([size(M0) 10]);
for ii=1:10
V(:,:,:,ii) = mrQ_PolyFitOrder(Mask,M0,ii);
end

for ii=1:10
    Coil=V(:,:,:,ii);
PErr(ii)=100*std( (M0(Mask)-Coil(Mask))./M0(Mask));
end

figure;
plot([1:10],PErr,'-k*');
grid on
set(gca,'YTick',[ 5 10 15 20 25 30 ])
    ylabel('Percent error','FontSize',20);
    xlabel('Polynomial order','FontSize',20);

