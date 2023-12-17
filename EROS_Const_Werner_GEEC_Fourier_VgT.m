clear;
fprintf('\n ###################################################### ');
fprintf('\n %s',datetime);
%% Model parameters, computation grid, density
%%%%%%%% EROS model
load EROS
nVpF=size(eros11272_22540,1);
nFace=22540;nVert=nVpF-nFace;
Vert=eros11272_22540(1:nVert,:);
Faces=eros11272_22540(nVert+1:end,:)+1; % Note the index start from 0 for EROS data
%%%%%%%% The external cuboid of the model
X1=min(Vert(:,1));X2=max(Vert(:,1));
Y1=min(Vert(:,2));Y2=max(Vert(:,2));
Z1=min(Vert(:,3));Z2=max(Vert(:,3));
%%%%%%%% Computation grid
% < ---------------------------- Parameter ---------------------------- > %
xmin=-20;dx=0.1;xmax=20;
ymin=-10;dy=0.1;ymax=10;
z0=-6.3;
% < ------------------------------ End  ------------------------------- > %
xgv=xmin:dx:xmax;ygv=ymin:dy:ymax;
nx=length(xgv);ny=length(ygv);
[X2d,Y2d]=meshgrid(xgv,ygv);
%%%%%%%% Density
% < ---------------------------- Parameter ---------------------------- > %
rho0=2670;
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Model parameters...');
fprintf('\n Number of vertices: %d.',nVert);
fprintf('\n Number of faces: %d.',nFace);
fprintf('\n Density: %.2f kg/m^3.',rho0);
fprintf('\n The external cuboid of the model: X:[%.2f,%.2f] km; Y:[%.2f,%.2f] km; Z:[%.2f,%.2f] km',...
    X1,X2,Y1,Y2,Z1,Z2);
fprintf('\n Computation plane altitude: %.2f km',z0);
fprintf('\n dx: %.2f km; dy: %.2f km; nx: %d; ny: %d; total computation points: %d',dx,dy,nx,ny,nx*ny);
%% Start computaion
%%%%%%%% Analytical solution of Saraswati_2019_JG
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Analytical solution of Saraswati_2019_JG...');
% Sample code to use GEEC %
% tic;
% Gcalc='grad';
% coord_calc=[X2d(:),Y2d(:),z0*ones(size(X2d(:)))];
% gTs_GEEC=geec(coord_calc,Vert,Faces,rho0,Gcalc);
% gTs_GEEC=6.673e-11/6.67408e-11*gTs_GEEC; % For different G
% gx_GEEC=1e3*reshape(gTs_GEEC(:,4),ny,nx);gy_GEEC=1e3*reshape(gTs_GEEC(:,5),ny,nx);gz_GEEC=1e3*reshape(gTs_GEEC(:,6),ny,nx);
% Txx_GEEC=reshape(gTs_GEEC(:,7),ny,nx); Txy_GEEC=reshape(gTs_GEEC(:,8),ny,nx); Txz_GEEC=reshape(gTs_GEEC(:,9),ny,nx);
% Tyy_GEEC=reshape(gTs_GEEC(:,10),ny,nx);Tyz_GEEC=reshape(gTs_GEEC(:,11),ny,nx);Tzz_GEEC=reshape(gTs_GEEC(:,12),ny,nx);
% time_gT_Saraswati=toc;
%%%% load result
matFileSaraswati=sprintf('EROS_Const%d_%d_%d_z0%g_dx%g_dy%g_gT_Saraswati',...
    rho0,nVert,nFace,z0,dx,dy);
matFileSaraswati = replace(matFileSaraswati,{'+','-','.'},{'p','m','d'});
load(matFileSaraswati);
%%%% Print the total number of computation points and time cost
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_gT_Saraswati);
fprintf('\n Number of facets computed per second: %d faces/sec',nx*ny*nFace/time_gT_Saraswati);
%%%%%%%% Werner_1997_CMDA
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Analytical solution of Werner_1997_CMDA...');
%%%% load result
matFileWerner=sprintf('EROS_Const%d_%d_%d_z0%g_dx%g_dy%g_VgT_Werner',...
    rho0,nVert,nFace,z0,dx,dy);
matFileWerner = replace(matFileWerner,{'+','-','.'},{'p','m','d'});
load(matFileWerner);
%%%% Print the total number of computation points and time cost
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_dVgT_Werner);
fprintf('\n Number of facets computed per second: %d faces/sec',nx*ny*nFace/time_dVgT_Werner);
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
t1=tic;
[dV_gft,gx_gft,gy_gft,gz_gft,Txx_gft,Tyy_gft,Tzz_gft,Txy_gft,Txz_gft,Tyz_gft]...
    =phi_xyz_PolyhedronSI_Parker(xgv,ygv,z0,...
    Vert,Faces,rho0,NT,Mx,My,Ntheta,Nk,eps_f,eps_b,R0,R1,typeNUFFT);
time_polyParker=toc(t1);
%%%% Print algorithm parameters and time cost
fprintf('\n NUFFT-Parker algorithm parameter:');
fprintf('\n kappa: %.2f; zmin: %.4f km; zmax: %.4f km; zref: %.4f km; chi: %.2f; N_T: %d',kappa,Z1,Z2,zref,chi,NT);
fprintf('\n (Mx=My): %d; Ntheta: %d; Nk: %d;eps: %g; R0: %d; R1: %d',...
    Mx,Ntheta,Nk,eps_f,R0,R1);
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_polyParker);
fprintf('\n Number of facets computed per second: %d faces/sec',nx*ny*nFace/time_polyParker);
%% Mapping
%%%%%%%% Werner_1997_CMDA as reference, to observe the error of Fourier method
del_gft_dV=dV_gft-dV_Werner;
del_gft_gx=gx_gft-gx_Werner;del_gft_gy=gy_gft-gy_Werner;del_gft_gz=gz_gft-gz_Werner;
del_gft_Txx=Txx_gft-Txx_Werner;del_gft_Tyy=Tyy_gft-Tyy_Werner;del_gft_Tzz=Tzz_gft-Tzz_Werner;
del_gft_Txy=Txy_gft-Txy_Werner;del_gft_Txz=Txz_gft-Txz_Werner;del_gft_Tyz=Tyz_gft-Tyz_Werner;
Datas_Werner_Fourier={dV_Werner,dV_gft,del_gft_dV;
    gx_Werner,gx_gft,del_gft_gx;
    gy_Werner,gy_gft,del_gft_gy;
    gz_Werner,gz_gft,del_gft_gz;
    Txx_Werner,Txx_gft,del_gft_Txx;
    Tyy_Werner,Tyy_gft,del_gft_Tyy;
    Tzz_Werner,Tzz_gft,del_gft_Tzz;
    Txy_Werner,Txy_gft,del_gft_Txy;
    Txz_Werner,Txz_gft,del_gft_Txz;
    Tyz_Werner,Tyz_gft,del_gft_Tyz
    };
%%%%%%%% Werner_1997_CMDA as reference, to observe the error of GEEC method
del_GEEC_gx=gx_GEEC-gx_Werner;del_GEEC_gy=gy_GEEC-gy_Werner;del_GEEC_gz=gz_GEEC-gz_Werner;
del_GEEC_Txx=Txx_GEEC-Txx_Werner;del_GEEC_Tyy=Tyy_GEEC-Tyy_Werner;del_GEEC_Tzz=Tzz_GEEC-Tzz_Werner;
del_GEEC_Txy=Txy_GEEC-Txy_Werner;del_GEEC_Txz=Txz_GEEC-Txz_Werner;del_GEEC_Tyz=Tyz_GEEC-Tyz_Werner;
Datas_Werner_GEEC={[],[],[]; % for convenience
    gx_Werner,gx_GEEC,del_GEEC_gx;
    gy_Werner,gy_GEEC,del_GEEC_gy;
    gz_Werner,gz_GEEC,del_GEEC_gz;
    Txx_Werner,Txx_GEEC,del_GEEC_Txx;
    Tyy_Werner,Tyy_GEEC,del_GEEC_Tyy;
    Tzz_Werner,Tzz_GEEC,del_GEEC_Tzz;
    Txy_Werner,Txy_GEEC,del_GEEC_Txy;
    Txz_Werner,Txz_GEEC,del_GEEC_Txz;
    Tyz_Werner,Tyz_GEEC,del_GEEC_Tyz
    };
% < ---------------------------- Parameter ---------------------------- > %
fs1=10;fs2=12;
showFlag=1; % 1: Show Fourier-Werner; 2: Show GEEC-Werner
% < ------------------------------ End  ------------------------------- > %
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
%%%%%%%% Potential and vector gravity fields
figure1= figure('Color',[1 1 1]);
if(showFlag==1)
    pickGrav=[1,2,3,4];
    Datas=Datas_Werner_Fourier;
else
    pickGrav=[2,3,4];
    Datas=Datas_Werner_GEEC;
end
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
%%%%%%%% Difference between Fourier and Werner_1997_CMDA solutions
stas_gft_dV=sta_min_max_mean_rms(dV_gft,dV_Werner,cutEPS_dVg);
stas_gft_gx=sta_min_max_mean_rms(gx_gft,gx_Werner,cutEPS_dVg);
stas_gft_gy=sta_min_max_mean_rms(gy_gft,gy_Werner,cutEPS_dVg);
stas_gft_gz=sta_min_max_mean_rms(gz_gft,gz_Werner,cutEPS_dVg);
stas_gft_Txx=sta_min_max_mean_rms(Txx_gft,Txx_Werner,cutEPS_T);
stas_gft_Tyy=sta_min_max_mean_rms(Tyy_gft,Tyy_Werner,cutEPS_T);
stas_gft_Tzz=sta_min_max_mean_rms(Tzz_gft,Tzz_Werner,cutEPS_T);
stas_gft_Txy=sta_min_max_mean_rms(Txy_gft,Txy_Werner,cutEPS_T);
stas_gft_Txz=sta_min_max_mean_rms(Txz_gft,Txz_Werner,cutEPS_T);
stas_gft_Tyz=sta_min_max_mean_rms(Tyz_gft,Tyz_Werner,cutEPS_T);
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Difference between Fourier and Werner_1997_CMDA solutions: ');
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


