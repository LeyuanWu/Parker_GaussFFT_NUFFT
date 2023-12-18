function [dV,gx,gy,gz,Txx,Tyy,Tzz,Txy,Txz,Tyz]=...
    phi_xyz_PolyhedronVI_Parker(xgv,ygv,z0,Nodes,Elements,rho,...
    NT,Mx,My,Ntheta,Nk,eps_f,eps_b,R0,R1,typeNUFFT)
% ----------------------------------------------------------------------- %
% M-file to calculate the complete gravity caused by a polyhedron with
% variable density using Modified Parker's method
% ----------------------------------------------------------------------- %
% ---- input ---- %
% xgv,ygv,z0: The 2D computation plane
% Nodes,Elements: the tetrahedron mesh
% rho: 3D density function
% Mx,My: Over-sampling of the Gauss-FFT algorithm
% Ntheta,Nk: Sampling of the spherical NUFFT algorithm
% eps_f,eps_b: precision request of the forward and backward NUFFT algorithm
% R0,R1: parameters for the partition function p_2(k)
% tyeNUFFT: string, 'ciamnufft' for codes at https://cims.nyu.edu/cmcl/nufft/nufft.html
%                   'finufft' for faster codes at https://finufft.readthedocs.io/en/latest/
% ---- output ---- %
% dV: gravity potential
% gx,gy,gz: gravity vector
% Txx,Tyy,Tzz,Txy,Txz,Tyz: gravity gradient tensor
% ---- Units ----%
% length in km, density in kg/m^3,
% gravity potential in m^2/s^2, vector in mGal, tensor in 1e-9/s^2 (Eotvos)
% ----------------------------------------------------------------------- %
%%%%%%%% Constants
G=6.673e-11;
CdV=2*pi*G*1e6;
Cg=2*pi*G*1e5*1e3;
Cten=2*pi*G*1e9;
%%%%%%%% Shift XY coordinates to center at the origin (x,y)=(0,0)
meanX=mean(xgv);xgv=xgv-meanX;Nodes(:,1)=Nodes(:,1)-meanX;
meanY=mean(ygv);ygv=ygv-meanY;Nodes(:,2)=Nodes(:,2)-meanY;
%%%%%%%% Grid parameters
nx=length(xgv);ny=length(ygv);
dx=xgv(2)-xgv(1);dy=ygv(2)-ygv(1);
nx1=ceil(-nx/2);nx2=ceil(nx/2-1);
ny1=ceil(-ny/2);ny2=ceil(ny/2-1);
NX=(nx1:1:nx2)';NY=(ny1:1:ny2)';
dkx=2*pi/(nx*dx);dky=2*pi/(ny*dy);
kxgv0=dkx*NX;kygv0=dky*NY;
[KX02d,KY02d]=meshgrid(kxgv0,kygv0);
alpha=xgv(1)/dx-nx1;beta=ygv(1)/dy-ny1;
%%%%%%%% VI over each tetrahedral element
fprintf('\n ******************** ');
fprintf('\n Construct the mass point model...');
tf1=tic;
Ras=Nodes(Elements(:,1),:);Rbs=Nodes(Elements(:,2),:);
Rcs=Nodes(Elements(:,3),:);Rds=Nodes(Elements(:,4),:);
clear Nodes Elements;
Qs=(Ras+Rbs+Rcs+Rds)/4;
Mabs=sqrt(sum((Rbs-Ras).^2,2));Macs=sqrt(sum((Rcs-Ras).^2,2));Mbcs=sqrt(sum((Rcs-Rbs).^2,2));
Mads=sqrt(sum((Rds-Ras).^2,2));Mbds=sqrt(sum((Rds-Rbs).^2,2));Mcds=sqrt(sum((Rds-Rcs).^2,2));
Ds=max([Mabs,Macs,Mbcs,Mads,Mbds,Mcds],[],2);
clear Mabs Macs Mbcs Mads Mbds Mcds;
Taus=abs(Qs(:,3)-z0)./Ds; % A simple way to get R when z0 is constant
threTaus=[inf,20,2,1,0];
nTaus=length(threTaus)-1;
cX=cell(nTaus,1);cY=cell(nTaus,1);cZ=cell(nTaus,1);
cW=cell(nTaus,1);
for iTau=1:1:nTaus
    pInd=((Taus>=threTaus(iTau+1)) & (Taus<threTaus(iTau)));
    if(any(pInd))
        pRas=Ras(pInd,:);pRbs=Rbs(pInd,:);pRcs=Rcs(pInd,:);pRds=Rds(pInd,:);
        N=iTau;
        [tXs,tYs,tZs,tWs]=tetraquad(pRas,pRbs,pRcs,pRds,N);
        cX{iTau}=tXs(:);cY{iTau}=tYs(:);cZ{iTau}=tZs(:);cW{iTau}=tWs(:);
    else
        continue;
    end
end
Xsrc=cell2mat(cX);Ysrc=cell2mat(cY);Zsrc=cell2mat(cZ);Wsrc=cell2mat(cW);
clear tXs tYs tZs tWs cX cY cZ cW;
Rho=rho(Xsrc,Ysrc,Zsrc);
WsrcRho=Wsrc.*Rho;
clear Wsrc Rho;
nsrc=length(Xsrc);
time_mp=toc(tf1);
fprintf('\n Total mass points: %d; Time cost: %.2f sec',nsrc,time_mp);
%%%%%%%% Optimal reference Zs levels for optimal convergence
minDepth=min(Zsrc)-z0;maxDepth=max(Zsrc)-z0;
nLevel=ceil(log(maxDepth/minDepth)/log(3));
Levels=(0:1:nLevel-1)';
SplitZs=z0+[minDepth*3.^Levels;1.001*maxDepth];
%%%%%%%% Split the mass points into several layers
XsrcLs=cell(nLevel,1);YsrcLs=cell(nLevel,1);ZsrcLs=cell(nLevel,1);WRhos=cell(nLevel,1);
for iL=1:1:nLevel
    pInd=Zsrc>=SplitZs(iL) & Zsrc<SplitZs(iL+1);
    XsrcLs{iL}=Xsrc(pInd);YsrcLs{iL}=Ysrc(pInd);ZsrcLs{iL}=Zsrc(pInd);WRhos{iL}=WsrcRho(pInd);
end
clear Xsrc Ysrc Zsrc WsrcRho;
%%%%%%%% Prepare R0 and R1
R0=R0*min(dkx,dky);R1=R1*min(dkx,dky); % From 2014Jiang_SIAM
%%%%%%%% Gauss-FFT to compute the regular part I1
dV_I1=zeros(ny,nx);
gx_I1=zeros(ny,nx);gy_I1=zeros(ny,nx);gz_I1=zeros(ny,nx);
Txx_I1=zeros(ny,nx);Tyy_I1=zeros(ny,nx);Tzz_I1=zeros(ny,nx);
Txy_I1=zeros(ny,nx);Txz_I1=zeros(ny,nx);Tyz_I1=zeros(ny,nx);
fprintf('\n ******************** ');
fprintf('\n Start calculating the Gauss-FFT part...');
tf2=tic;
%%%% Gauss weights and nodes
[Ax,Wx]=lgwt(Mx,0,1);[Ay,Wy]=lgwt(My,0,1);
if(mod(My,2)==0)
    n2=length(Wy)/2;Wy=2*Wy(1:n2);Ay=Ay(1:n2);
else
    n2=(length(Wy)+1)/2;Wy=Wy(1:n2);Wy(1:n2-1)=2*Wy(1:n2-1);Ay=Ay(1:n2);
end
[AX,AY]=meshgrid(Ax,Ay);AX=AX(:);AY=AY(:);
[WX,WY]=meshgrid(Wx,Wy);WX=WX(:);WY=WY(:);
%%%% Gauss-FFT summation
nGFT=length(AX);
for iGFT=1:1:nGFT
    ksi=AX(iGFT);eta=AY(iGFT);
    KX2d=KX02d+ksi*dkx;KY2d=KY02d+eta*dky;
    K2d=sqrt(KX2d.^2+KY2d.^2);
    KConstI1=exp(z0*K2d)./K2d;
    ThetaI1_v=zeros(ny,nx);
    for iL=1:1:nLevel
        delta=(SplitZs(iL)+SplitZs(iL+1))/2;
        tXsrc=XsrcLs{iL};tYsrc=YsrcLs{iL};tZsrc=ZsrcLs{iL};
        tWRho=WRhos{iL};ntsrc=numel(tXsrc);
        shift_src=exp(-1i*ksi*dkx*tXsrc).*exp(-1i*eta*dky*tYsrc);
        temp=zeros(ny,nx);
        theta_v_n=tWRho.*shift_src;
        for n=0:1:NT
            if(strcmp(typeNUFFT,'ciamnufft'))
                Theta_v_n=ntsrc*nufft2d1(ntsrc,dky*tYsrc,dkx*tXsrc,theta_v_n,-1,eps_f,ny,nx);
            else
                fi_opt.debug = 0; fi_opt.upsampfac=0;fi_opt.chkbnds=0;
                FFTW_ESTIMATE = bitshift(1,6); fi_opt.fftw = FFTW_ESTIMATE;
                Theta_v_n=finufft2d1(dky*tYsrc,dkx*tXsrc,theta_v_n,-1,eps_f/10,ny,nx,fi_opt);
            end
            temp=temp+(-K2d).^n/factorial(n).*Theta_v_n;
            theta_v_n=theta_v_n.*(tZsrc-delta);
        end
        ThetaI1_v=ThetaI1_v+exp(-delta*K2d).*temp;
    end
    P2k2d=1/2*(1-erf(12*(K2d-(R0+R1)/2)./(R1-R0)));
    ThetaI1_v=(1-P2k2d).*KConstI1.*ThetaI1_v;
    ThetaI1_x=1i*KX2d.*ThetaI1_v;
    ThetaI1_y=1i*KY2d.*ThetaI1_v;
    ThetaI1_z=K2d.*ThetaI1_v;
    ThetaI1_xx=-KX2d.^2.*ThetaI1_v;
    ThetaI1_yy=-KY2d.^2.*ThetaI1_v;
    ThetaI1_zz=K2d.^2.*ThetaI1_v;
    ThetaI1_xy=-KX2d.*KY2d.*ThetaI1_v;
    ThetaI1_xz=1i*KX2d.*K2d.*ThetaI1_v;
    ThetaI1_yz=1i*KY2d.*K2d.*ThetaI1_v;
    dV_I1 =dV_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_v,ksi,eta,alpha,beta));
    gx_I1 =gx_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_x,ksi,eta,alpha,beta));
    gy_I1 =gy_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_y,ksi,eta,alpha,beta));
    gz_I1 =gz_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_z,ksi,eta,alpha,beta));
    Txx_I1 =Txx_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_xx,ksi,eta,alpha,beta));
    Tyy_I1 =Tyy_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_yy,ksi,eta,alpha,beta));
    Tzz_I1 =Tzz_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_zz,ksi,eta,alpha,beta));
    Txy_I1 =Txy_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_xy,ksi,eta,alpha,beta));
    Txz_I1 =Txz_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_xz,ksi,eta,alpha,beta));
    Tyz_I1 =Tyz_I1+WX(iGFT)*WY(iGFT)...
        *real(sifft2(ThetaI1_yz,ksi,eta,alpha,beta));
end
dV_I1=CdV/(dx*dy)*dV_I1;
tCg=Cg/(dx*dy);
gx_I1=tCg*gx_I1;gy_I1=tCg*gy_I1;gz_I1=tCg*gz_I1;
tCten=Cten/(dx*dy);
Txx_I1=tCten*Txx_I1;Tyy_I1=tCten*Tyy_I1;Tzz_I1=tCten*Tzz_I1;
Txy_I1=tCten*Txy_I1;Txz_I1=tCten*Txz_I1;Tyz_I1=tCten*Tyz_I1;
time_GFT=toc(tf2);
fprintf('\n Time cost: %.2f sec',time_GFT);
%%%%%%%% NUFFT to compute the singular part I2
fprintf('\n ******************** ');
fprintf('\n Start calculating the NUFFT part...');
tf3=tic;
%%%% Polar NUFFT sampling
thetagv=linspace(0,2*pi-2*pi/Ntheta,Ntheta);
wtheta=2*pi/Ntheta*ones(Ntheta,1);
[magkgv,wk]=lgwt(Nk,0,R1);
[THETA,MAGK]=meshgrid(thetagv,magkgv);
THETA=THETA(:);MAGK=MAGK(:);
[WTHETA,WK]=meshgrid(wtheta,wk);
Ws=WTHETA.*WK;Ws=Ws(:);
KX_NUFFT=MAGK.*cos(THETA);
KY_NUFFT=MAGK.*sin(THETA);
K_NUFFT=sqrt(KX_NUFFT.^2+KY_NUFFT.^2);
%%%% Number of sampling and shift parameters
nspec=Ntheta*Nk;
shift_spec=exp(1i*alpha*dx*KX_NUFFT).*exp(1i*beta*dy*KY_NUFFT);
%%%% The Non-uniform spectrum using NUFFT type 3
ThetaI2_v=zeros(nspec,1);
for iL=1:1:nLevel
    delta=(SplitZs(iL)+SplitZs(iL+1))/2;
    tXsrc=XsrcLs{iL};tYsrc=YsrcLs{iL};tZsrc=ZsrcLs{iL};
    tWRho=WRhos{iL};ntsrc=numel(tXsrc);
    temp=zeros(nspec,1);
    for n=0:1:NT
        if(strcmp(typeNUFFT,'ciamnufft'))
            Theta_v_n=nufft2d3(ntsrc,dy*tYsrc,dx*tXsrc,tWRho,-1,eps_f,nspec,KY_NUFFT/dy,KX_NUFFT/dx);
        else
            Theta_v_n=finufft2d3(dy*tYsrc,dx*tXsrc,tWRho,-1,eps_f/10,KY_NUFFT/dy,KX_NUFFT/dx,fi_opt);
        end
        temp=temp+(-K_NUFFT).^n/factorial(n).*Theta_v_n;
        tWRho=tWRho.*(tZsrc-delta);
    end
    ThetaI2_v=ThetaI2_v+exp(-delta*K_NUFFT).*temp;
end
P2k_NUFFT=1/2*(1-erf(12*(K_NUFFT-(R0+R1)/2)./(R1-R0)));
ThetaI2_v=Ws.*P2k_NUFFT.*exp(z0*K_NUFFT).*ThetaI2_v;
ThetaI2_x=1i*KX_NUFFT.*ThetaI2_v;
ThetaI2_y=1i*KY_NUFFT.*ThetaI2_v;
ThetaI2_z=K_NUFFT.*ThetaI2_v;
ThetaI2_xx=-KX_NUFFT.^2.*ThetaI2_v;
ThetaI2_yy=-KY_NUFFT.^2.*ThetaI2_v;
ThetaI2_zz=K_NUFFT.^2.*ThetaI2_v;
ThetaI2_xy=-KX_NUFFT.*KY_NUFFT.*ThetaI2_v;
ThetaI2_xz=1i*KX_NUFFT.*K_NUFFT.*ThetaI2_v;
ThetaI2_yz=1i*KY_NUFFT.*K_NUFFT.*ThetaI2_v;
%%%% The gridded gravity field using NUFFT type 1
if(strcmp(typeNUFFT,'ciamnufft'))
    dV_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_v,1,eps_b,ny,nx);
    gx_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_x,1,eps_b,ny,nx);
    gy_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_y,1,eps_b,ny,nx);
    gz_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_z,1,eps_b,ny,nx);
    Txx_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xx,1,eps_b,ny,nx);
    Tyy_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_yy,1,eps_b,ny,nx);
    Tzz_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_zz,1,eps_b,ny,nx);
    Txy_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xy,1,eps_b,ny,nx);
    Txz_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xz,1,eps_b,ny,nx);
    Tyz_I2=nspec*nufft2d1(nspec,dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_yz,1,eps_b,ny,nx);
else
    dV_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_v,1,eps_b/10,ny,nx,fi_opt);
    gx_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_x,1,eps_b/10,ny,nx,fi_opt);
    gy_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_y,1,eps_b/10,ny,nx,fi_opt);
    gz_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_z,1,eps_b/10,ny,nx,fi_opt);
    Txx_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xx,1,eps_b/10,ny,nx,fi_opt);
    Tyy_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_yy,1,eps_b/10,ny,nx,fi_opt);
    Tzz_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_zz,1,eps_b/10,ny,nx,fi_opt);
    Txy_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xy,1,eps_b/10,ny,nx,fi_opt);
    Txz_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_xz,1,eps_b/10,ny,nx,fi_opt);
    Tyz_I2=finufft2d1(dy*KY_NUFFT,dx*KX_NUFFT,shift_spec.*ThetaI2_yz,1,eps_b/10,ny,nx,fi_opt);
end
dV_I2=CdV/(2*pi)^2*real(dV_I2);
tCg=Cg/(2*pi)^2;
gx_I2=tCg*real(gx_I2);gy_I2=tCg*real(gy_I2);gz_I2=tCg*real(gz_I2);
tCten=Cten/(2*pi)^2;
Txx_I2=tCten*real(Txx_I2);Tyy_I2=tCten*real(Tyy_I2);Tzz_I2=tCten*real(Tzz_I2);
Txy_I2=tCten*real(Txy_I2);Txz_I2=tCten*real(Txz_I2);Tyz_I2=tCten*real(Tyz_I2);
time_NUFFT=toc(tf3);
fprintf('\n Time cost: %.2f sec',time_NUFFT);
%%%%%%%% Adding the two parts together
dV=dV_I1+dV_I2;
gx=gx_I1+gx_I2;gy=gy_I1+gy_I2;gz=gz_I1+gz_I2;
Txx=Txx_I1+Txx_I2;Tyy=Tyy_I1+Tyy_I2;Tzz=Tzz_I1+Tzz_I2;
Txy=Txy_I1+Txy_I2;Txz=Txz_I1+Txz_I2;Tyz=Tyz_I1+Tyz_I2;
fprintf('\n Complete!');
