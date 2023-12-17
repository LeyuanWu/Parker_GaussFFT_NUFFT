function [dV,gx,gy,gz,Txx,Tyy,Tzz,Txy,Txz,Tyz]=...
    phi_xyz_sphere(xgv,ygv,z0,xc,yc,zc,a,rho)
% -------------------------------------------------------------------------
% M-file to calculate the complete gravity caused by a uniform sphere.
% Only valid for computation points outside the sphere.
% -------------------------------------------------------------------------
% ---- input ---- %
% xgv,ygv,z0: The 2D computation plane
% Note: if xgv,ygv,z0 are of the same length, they represent  irregular computation
%       points; if not, they represent a computation plane with z0 the altitude
% xc,yc,zc: Center of the sphere
% a: Radius of the sphere
% rho: density
% ---- output ---- %
% dV: gravity potential
% gx,gy,gz: gravity vector
% Txx,Tyy,Tzz,Txy,Txz,Tyz: gravity gradient tensor
% ---- Units ----%
% length in km, density in kg/m^3, 
% gravity potential in m^2/s^2, vector in mGal, tensor in 1e-9/s^2 (Eotvos)
% -------------------------------------------------------------------------
%%%%%%%% Determine if the computation is carried out on a horizontal plane
if(length(xgv)~=length(z0))
    flag=1;
    nx=length(xgv);ny=length(ygv);
    [X2d,Y2d]=meshgrid(xgv,ygv);
    Z2d=z0*ones(size(X2d));
    XPs=X2d(:);YPs=Y2d(:);ZPs=Z2d(:);
else
    flag=2;
    XPs=xgv(:);YPs=ygv(:);ZPs=z0(:);
end
%%%%%%%% Constants
gamma=6.673e-11;
si2mg=1.e5;km2m=1.e3;
tmass=4.*pi*rho*(a^3)/3;
CdV=gamma*tmass*km2m^2;
Cg=-gamma*tmass*si2mg*km2m;
Cten=-gamma*tmass*1e9;
%%%%%%%% Computation
nP=length(XPs);
dV=zeros(nP,1);
gx=zeros(nP,1);gy=zeros(nP,1);gz=zeros(nP,1);
Txx=zeros(nP,1);Tyy=zeros(nP,1);Tzz=zeros(nP,1);
Txy=zeros(nP,1);Txz=zeros(nP,1);Tyz=zeros(nP,1);
for iP=1:1:nP
    xp=XPs(iP);
    yp=YPs(iP);
    zp=ZPs(iP);
    rz=zp-zc;
    rx=xp-xc;
    ry=yp-yc;
    r=sqrt(rx^2+ry^2+rz^2);
    r2=r^2;
    r3=r^3;
    r5=r^5;
    dV(iP)=1/r;
    gx(iP)=rx/r3;
    gy(iP)=ry/r3;
    gz(iP)=rz/r3;
    Txx(iP)=(r2-3*rx^2)/r5;
    Tyy(iP)=(r2-3*ry^2)/r5;
    Tzz(iP)=(r2-3*rz^2)/r5;
    Txy(iP)=-3*(rx*ry)/r5;
    Txz(iP)=-3*(rx*rz)/r5;
    Tyz(iP)=-3*(ry*rz)/r5;
end
dV=CdV*dV; % m^2/s^2
gx=Cg*gx;gy=Cg*gy;gz=Cg*gz; % mGal
Txx=Cten*Txx;Tyy=Cten*Tyy;Tzz=Cten*Tzz; % Eotvos
Txy=Cten*Txy;Txz=Cten*Txz;Tyz=Cten*Tyz;
%%%%%%%% Transform to 2D plane anomalies
if(flag==1)
    dV=reshape(dV,ny,nx);
    gx=reshape(gx,ny,nx);
    gy=reshape(gy,ny,nx);
    gz=reshape(gz,ny,nx);
    Txx=reshape(Txx,ny,nx);
    Tyy=reshape(Tyy,ny,nx);
    Tzz=reshape(Tzz,ny,nx);
    Txy=reshape(Txy,ny,nx);
    Txz=reshape(Txz,ny,nx);
    Tyz=reshape(Tyz,ny,nx);
end




