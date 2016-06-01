function MappimngErr=mrQ_testCSFFit(B1file,rawDatfile,Gfile,PDfile,T1file)


load(rawDatfile);
B1=readFileNifti(B1file);B1=double(B1.data);
G=readFileNifti(Gfile);B1=double(G.data);
T1=readFileNifti(T1file);T1=double(T1.data);
PD=readFileNifti(PDfile);B1=double(PD.data);

tr=S.tr;
fa=S.fipAngle;

for ii=1:length(tr)
    fa=flipAngles(ii).*B1;
fa = fa./180.*pi;

Sc(:,ii) =G(:)(1-exp(-tr(ii)./T1(:))).*sin(fa(:))./(1-exp(-tr(ii)./T1(:)).*cos(fa));

end

Sig=


MappimngErr=Sig-Sc;
