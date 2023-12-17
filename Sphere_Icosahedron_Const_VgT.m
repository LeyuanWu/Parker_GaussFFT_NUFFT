clear;
fprintf('\n ###################################################### ');
fprintf('\n %s',datetime);
%% Spherical triangulation based on an Icosahedral
%%%%%%%% Icosahedral
r=1;ddV=(1+sqrt(5))/2;
V1=[0;0;0;0;-1;-1;1;1;-ddV;ddV;ddV;-ddV;];
V2=[-1;-1;1;1;-ddV;ddV;ddV;-ddV;0;0;0;0;];
V3=[-ddV;ddV;ddV;-ddV;0;0;0;0;-1;-1;1;1;];
FacesIco= [9 4 1;1 5 9;1 8 5;10 8 1;4 10 1;
    12 5 2; 12 2 3; 12 3 6; 12 6 9; 12 9 5;
    11 10 7; 11 8 10; 11 2 8; 11 3 2; 11 7 3;
    2 5 8; 10 4 7; 3 7 6; 6 7 4; 6 4 9;];
[THETA,PHI,~]=cart2sph(V1,V2,V3);
R=r.*ones(size(V1(:,1)));
[V1,V2,V3]=sph2cart(THETA,PHI,R);
VertIco=[V1 V2 V3];
%%%%%%%% Triangular subdivision
% < ---------------------------- Parameter ---------------------------- > %
Level=8;
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Spherical triangulation based on an Icosahedral: Level=%d...',Level);
tic;
N=2^Level;
nFace=size(FacesIco,1);
VertUnit=zeros(nFace*(N+1)*(N+2)/2,3);
Faces=zeros(nFace*N^2,3);
for iFace=1:1:nFace
    Rs=VertIco(FacesIco(iFace,:),:);
    Ra=Rs(1,:)';Rb=Rs(2,:)';Rc=Rs(3,:)';
    Mab=Rb-Ra;Mac=Rc-Ra;Mbc=Rc-Rb;
    ab=norm(Mab);ac=norm(Mac);bc=norm(Mbc);
    ni_temp=[Mab(2)*Mac(3)-Mab(3)*Mac(2);Mab(3)*Mac(1)-Mab(1)*Mac(3);Mab(1)*Mac(2)-Mab(2)*Mac(1)];
    ni=ni_temp/norm(ni_temp);
    indVert0=(iFace-1)*(N+1)*(N+2)/2;
    indVert=indVert0+1;
    for irow=1:1:N+1
        for icol=1:1:irow
            VertUnit(indVert,:)=Ra+(irow-1)/N*Mab+(icol-1)/N*Mbc;
            indVert=indVert+1;
        end
    end
    indFace0=(iFace-1)*N^2;
    for irow=1:1:N
        for icol=1:1:irow
            if(icol<irow)
                ind1=indVert0+irow*(irow-1)/2+icol;
                ind2=ind1+irow;ind3=ind2+1;ind4=ind1+1;
                indFace=indFace0+(irow-1)^2+2*(icol-1)+1;
                Faces(indFace,:)=[ind1,ind2,ind3];
                Faces(indFace+1,:)=[ind1,ind3,ind4];
            else
                ind1=indVert0+irow*(irow+1)/2;
                ind2=ind1+irow;ind3=ind2+1;
                indFace=indFace0+(irow-1)^2+2*(icol-1)+1;
                Faces(indFace,:)=[ind1,ind2,ind3];
            end
        end
    end
end
VL1=VertUnit(:,1);VL2=VertUnit(:,2);VL3=VertUnit(:,3);
time_PGIM=toc;
fprintf('\n Complete the subdivision.');
fprintf('\n Number of vertices: %d.',size(VertUnit,1));
fprintf('\n Number of faces: %d.',size(Faces,1));
fprintf('\n Time cost: %.2f sec',time_PGIM);
%% (Sphere <-> Icosahedron) Model parameters, computation grid
%%%%%%%% Computation grid
% < ---------------------------- Parameter ---------------------------- > %
xmin=-10;dx=0.1;xmax=10;
ymin=-10;dy=0.1;ymax=10;
z0=0;
% < ------------------------------ End  ------------------------------- > %
xgv=xmin:dx:xmax;ygv=ymin:dy:ymax;
nx=length(xgv);ny=length(ygv);
[X2d,Y2d]=meshgrid(xgv,ygv);
%%%%%%%% Position, depth and density of the sphere
% < ---------------------------- Parameter ---------------------------- > %
Radius=1.5;
xc=0;yc=0;zc=2;
rho0=1000;
% < ------------------------------ End  ------------------------------- > %
%%%%%%% Geometry of the Sphere <-> Icosahedron Model
[THETAL,PHIL,~]=cart2sph(VL1,VL2,VL3);
[VL1,VL2,VL3]=sph2cart(THETAL,PHIL,Radius);
Vert=[VL1,VL2,VL3];
Vert(:,1)=Vert(:,1)+xc;
Vert(:,2)=Vert(:,2)+yc;
Vert(:,3)=Vert(:,3)+zc;
%%%%%%%% The external cuboid of the model
X1=min(Vert(:,1));X2=max(Vert(:,1));
Y1=min(Vert(:,2));Y2=max(Vert(:,2));
Z1=min(Vert(:,3));Z2=max(Vert(:,3));
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Model parameters...');
fprintf('\n Density: %.2f kg/m^3.',rho0);
fprintf('\n The external cuboid of the model: X:[%.2f,%.2f] km; Y:[%.2f,%.2f] km; Z:[%.2f,%.2f] km',...
    X1,X2,Y1,Y2,Z1,Z2);
fprintf('\n Computation plane altitude: %.2f km',z0);
fprintf('\n dx: %.2f km; dy: %.2f km; nx: %d; ny: %d; total computation points: %d',dx,dy,nx,ny,nx*ny);
%% Start computaion
%%%%%%%% Analytical solution of a uniform sphere
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Analytical solution of a uniform sphere...');
t1=tic;
[dV,gx,gy,gz,Txx,Tyy,Tzz,Txy,Txz,Tyz]...
    =phi_xyz_sphere(xgv,ygv,z0,xc,yc,zc,Radius,rho0);
time_sphere=toc(t1);
%%%% Print the total number of computation points and time cost
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_sphere);
%%%%%%%% NUFFT-based Fourier-Parker method
%%%% Parameters related to the convergence behavior of Parker's formula
kappa=(Z1-z0)/max([dx,dy]);
zref=(Z1+Z2)/2;
chi=(zref-Z1)/(zref-z0);
%%%% Algorithmic parameters
% < ---------------------------- Parameter ---------------------------- > %
Mx=2;My=2;
Ntheta=floor(max(nx,ny));Nk=floor(max(nx,ny)/3);
R0=0;R1=floor(max(nx,ny)/8);
eps_f=1e-4;eps_b=1e-4;
NT=10;
typeNUFFT='finufft';   
% ciamnufft for codes at https://cims.nyu.edu/cmcl/nufft/nufft.html
% finufft for codes at https://finufft.readthedocs.io/en/latest/
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n NUFFT-based Fourier-Parker method using *%s*...',typeNUFFT);
t2=tic;
[dV_gft,gx_gft,gy_gft,gz_gft,Txx_gft,Tyy_gft,Tzz_gft,Txy_gft,Txz_gft,Tyz_gft]...
    =phi_xyz_PolyhedronSI_Parker(xgv,ygv,z0,...
    Vert,Faces,rho0,NT,Mx,My,Ntheta,Nk,eps_f,eps_b,R0,R1,typeNUFFT);
time_polyParker=toc(t2);
%%%% Print algorithm parameters and time cost
fprintf('\n NUFFT-Parker algorithm parameter:');
fprintf('\n kappa: %.2f; zmin: %.4f km; zmax: %.4f km; zref: %.4f km; chi: %.2f; N_T: %d',kappa,Z1,Z2,zref,chi,NT);
fprintf('\n (Mx=My): %d; Ntheta: %d; Nk: %d;eps: %g; R0: %d; R1: %d',...
    Mx,Ntheta,Nk,eps_f,R0,R1);
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_polyParker);
%% Mapping
%%%%%%%% Difference between Fourier and Analytical solution
del_gft_dV=dV_gft-dV;
del_gft_gx=gx_gft-gx;del_gft_gy=gy_gft-gy;del_gft_gz=gz_gft-gz;
del_gft_Txx=Txx_gft-Txx;del_gft_Tyy=Tyy_gft-Tyy;del_gft_Tzz=Tzz_gft-Tzz;
del_gft_Txy=Txy_gft-Txy;del_gft_Txz=Txz_gft-Txz;del_gft_Tyz=Tyz_gft-Tyz;
Datas={dV,dV_gft,del_gft_dV;
    gx,gx_gft,del_gft_gx;
    gy,gy_gft,del_gft_gy;
    gz,gz_gft,del_gft_gz;
    Txx,Txx_gft,del_gft_Txx;
    Tyy,Tyy_gft,del_gft_Tyy;
    Tzz,Tzz_gft,del_gft_Tzz;
    Txy,Txy_gft,del_gft_Txy;
    Txz,Txz_gft,del_gft_Txz;
    Tyz,Tyz_gft,del_gft_Tyz
    };
Titles={'$V$','$\hat{V}$','$\delta(V)$';
    '$g_x$','$\hat{g_x}$','$\delta(g_x)$';
    '$g_y$','$\hat{g_y}$','$\delta(g_y)$';
    '$g_z$','$\hat{g_z}$','$\delta(g_z)$';
    '$T_{xx}$','$\hat{T_{xx}}$','$\delta(T_{xx})$';
    '$T_{yy}$','$\hat{T_{yy}}$','$\delta(T_{yy})$';
    '$T_{zz}$','$\hat{T_{zz}}$','$\delta(T_{zz})$';
    '$T_{xy}$','$\hat{T_{xy}}$','$\delta(T_{xy})$';
    '$T_{xz}$','$\hat{T_{xz}}$','$\delta(T_{xz})$';
    '$T_{yz}$','$\hat{T_{yz}}$','$\delta(T_{yz})$';
    };
% < ---------------------------- Parameter ---------------------------- > %
fs1=10;fs2=12;
% < ------------------------------ End  ------------------------------- > %
%%%%%%%% Potential and vector gravity fields
figure1= figure('Color',[1 1 1]);
pickGrav=[1,2,3,4];
nGrav=length(pickGrav);
for i=1:1:nGrav
    for j=1:1:3
        curGrav=pickGrav(i);
        subplot1=subplot(nGrav,3,3*(i-1)+j,'Parent',figure1,...
            'DataAspectRatio',[1 1 1],'FontSize',fs1);
        hold(subplot1,'all');grid(subplot1,'on');
        [c,h]=contourf(X2d,Y2d,Datas{curGrav,j});
        if(j==1)
            lvlStep=h.LevelStep;
        end
        if(j==2)
            h.LevelStep=lvlStep;
        end
        colormap(jet);colorbar;
        xlabel('x (km)');ylabel('y (km)');
        title(Titles{curGrav,j},'Interpreter','latex','FontSize',fs2);
    end
end
%%%%%%%% Gravity Gradient Tensor
figure2= figure('Color',[1 1 1]);
pickGrav=[5,6,7,8,9,10];
nGrav=length(pickGrav);
for i=1:1:nGrav
    for j=1:1:3
        curGrav=pickGrav(i);
        subplot1=subplot(nGrav,3,3*(i-1)+j,'Parent',figure2,...
            'DataAspectRatio',[1 1 1],'FontSize',fs1);
        hold(subplot1,'all');grid(subplot1,'on');
        [c,h]=contourf(X2d,Y2d,Datas{curGrav,j});
        if(j==1)
            lvlStep=h.LevelStep;
        end
        if(j==2)
            h.LevelStep=lvlStep;
        end
        colormap(jet);colorbar;
        xlabel('x (km)');ylabel('y (km)');
        title(Titles{curGrav,j},'Interpreter','latex','FontSize',fs2);
    end
end
%% Print information on screen (also for table)
% < ---------------------------- Parameter ---------------------------- > %
cutEPS_dVg=1e-3; % if value < cutEPS_dVg*max, set eps as nan
cutEPS_T=1e-3; % if value < cutEPS_T*max, set eps as nan
pickErrInfor=[1:8,13];
% < ------------------------------ End  ------------------------------- > %
%%%%%%%% Difference between Fourier and Sphere analytical solution
stas_gft_dV=sta_min_max_mean_rms(dV_gft,dV,cutEPS_dVg);
stas_gft_gx=sta_min_max_mean_rms(gx_gft,gx,cutEPS_dVg);
stas_gft_gy=sta_min_max_mean_rms(gy_gft,gy,cutEPS_dVg);
stas_gft_gz=sta_min_max_mean_rms(gz_gft,gz,cutEPS_dVg);
stas_gft_Txx=sta_min_max_mean_rms(Txx_gft,Txx,cutEPS_T);
stas_gft_Tyy=sta_min_max_mean_rms(Tyy_gft,Tyy,cutEPS_T);
stas_gft_Tzz=sta_min_max_mean_rms(Tzz_gft,Tzz,cutEPS_T);
stas_gft_Txy=sta_min_max_mean_rms(Txy_gft,Txy,cutEPS_T);
stas_gft_Txz=sta_min_max_mean_rms(Txz_gft,Txz,cutEPS_T);
stas_gft_Tyz=sta_min_max_mean_rms(Tyz_gft,Tyz,cutEPS_T);
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Difference between Fourier and Sphere analytical solution: ');
fprintf('\n min(ref) & max(ref) & mean(ref) & rms(ref) & min(err) & max(err) & mean(err) & rms(err) & E2')
fprintf('\n $V$ ');fprintf('& %.2e ',stas_gft_dV(pickErrInfor));fprintf(' \\\\');
fprintf('\n $g_x$ ');fprintf('& %.2e ',stas_gft_gx(pickErrInfor));fprintf(' \\\\');
fprintf('\n $g_y$ ');fprintf('& %.2e ',stas_gft_gy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $g_z$ ');fprintf('& %.2e ',stas_gft_gz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xx}$ ');fprintf('& %.2e ',stas_gft_Txx(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{yy}$ ');fprintf('& %.2e ',stas_gft_Tyy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{zz}$ ');fprintf('& %.2e ',stas_gft_Tzz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xy}$ ');fprintf('& %.2e ',stas_gft_Txy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xz}$ ');fprintf('& %.2e ',stas_gft_Txz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{yz}$ ');fprintf('& %.2e ',stas_gft_Tyz(pickErrInfor));fprintf(' \\\\');












