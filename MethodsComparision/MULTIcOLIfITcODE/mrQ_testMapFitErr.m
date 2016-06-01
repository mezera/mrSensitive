function MappimngErr=mrQ_testMapFitErr(B1file,rawDatfile,Gfile,PDfile,T1file,BMfile,boxsize)


load(rawDatfile);
B1=readFileNifti(B1file);B1=double(B1.data);
G=readFileNifti(Gfile);G=double(G.data);
T1=readFileNifti(T1file);T1=double(T1.data).*1000;
PD=readFileNifti(PDfile);PD=double(PD.data);
BM=readFileNifti(BMfile);BM=logical(BM.data);

flipAngles = [s(:).flipAngle];
    tr         = [s(:).TR];

for ii=1:length(tr)
    fa=flipAngles(ii).*B1;
fa = fa./180.*pi;


Sc(:,:,:,ii) =G.*PD.*(1-exp(-tr(ii)./T1)).*sin(fa)./(1-exp(-tr(ii)./T1).*cos(fa));

end
clear T1 B1 G PD

Sig=cat(4,s(:).imData);
clear s

MappimngErr=Sig-Sc;
MappimngErr1=Sig./Sc;

 
if notDefined('boxsize')
    boxsize(1)=30;
    boxsize(2)=40;
    boxsize(3)=20;
end

sz=size(Sig);sz=sz(1:3);
CSF1=ones(sz);
szH=round(sz./2);
XX=boxsize(1)./round(mmPerVox(1));
YY=boxsize(2)./round(mmPerVox(2));
ZZ=boxsize(3)./round(mmPerVox(3));

CSF1(szH(1)+XX:end,:,:)=0;
CSF1(1:szH(1)-XX,:,:)=0;

CSF1(:,1:szH(2)-YY,:)=0;
CSF1(:,szH(2)+YY:end,:,:)=0;

CSF1(:,:,1:szH(3)-ZZ)=0;
CSF1(:,:,szH(3)+ZZ:end)=0;


for ii=1:length(tr)
    dat=Sig(:,:,:,ii);
        fit=Sc(:,:,:,ii);

BMT=BM & dat>30;
BME= CSF1==1 & dat<30;
RateT(ii)=median(fit(BMT)./dat(BMT));
RateE(ii)=median(fit(BME)./dat(BME));


end

for ii=1:length(tr)
    Tmp= Sig(:,:,:,ii);
out=Sc(:,:,:,ii)==0 & Sig(:,:,:,ii)>0;
NoiseL(ii)=median(Tmp(out));
end

[tt,M01]    = relaxFitT1(cat(4,s(:).imData),flipAngles,tr(1),B1);



for ii=1:length(tr)

figure;

 subplot(2,2,1);imagesc(Sig(:,:,90,ii)); axis off; title(['signal flipangle = ' num2str(flipAngles(ii))  ]);colorbar 
  subplot(2,2,2);imagesc(Sc(:,:,90,ii)); axis off; title('fit');colorbar 
  subplot(2,2,3);imagesc(MappimngErr(:,:,90,ii)); axis off; title('sig-fit'); caxis([-5.5 5.50000]);colorbar 
  subplot(2,2,4);imagesc(MappimngErr1(:,:,90,ii)); axis off; title('sig/fit');caxis([.9 1.10000]);colorbar ; colormap hot
end




fa=flipAngles(1).*B1;
fa = fa./180.*pi;

T1csf=4300;
Im =G.*(1-exp(-tr(ii)./T1csf)).*sin(fa)./(1-exp(-tr(ii)./T1csf).*cos(fa));
PDcsf=Sig(:,:,:,1)./Im;

CSF2=CSF1 & T1>4000;

PDval=median((PDcsf(CSF2)));


err=(PDval.*Im-Sig(:,:,:,1))./Sig(:,:,:,1);

CSF3=abs(err)< std(err(CSF2)) & CSF2;
