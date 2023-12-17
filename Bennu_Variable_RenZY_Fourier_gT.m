clear;
fprintf('\n ###################################################### ');
fprintf('\n %s',datetime);
%% Bennu model decomposed into a tetrahedral mesh using TetGen
%%%%%%%% Bennu surface model
load Bennu
nVpF=size(Bennu25350_49152,1);
nFace=49152;nVert=nVpF-nFace;
Vert=Bennu25350_49152(1:nVert,:);
Faces=Bennu25350_49152(nVert+1:end,:);
%%%%%%%% The external cuboid of the model
X1=min(Vert(:,1));X2=max(Vert(:,1));
Y1=min(Vert(:,2));Y2=max(Vert(:,2));
Z1=min(Vert(:,3));Z2=max(Vert(:,3));
%%%%%%%% Bennu Tetrahedral mesh generated using the Bennu25350_49152 model
fileNode='Bennus49152_pq1.2_295272.node'; % 295272 Nodes
fileEle='Bennus49152_pq1.2_1432834.ele';  % 1432834 Elements
%%%%%%%% Processing tetrahedral mesh data
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Processing tetrahedral mesh data... ');
t1=tic;
fileID = fopen(fileNode);
formatSpec1 = '%d %d %d %d';sizeFirstLine = [4 1];
FirstLine= fscanf(fileID,formatSpec1,sizeFirstLine);
formatSpec2 = '%d %f %f %f';nNode=FirstLine(1);sizeNodes=[4,nNode];
Nodes = fscanf(fileID,formatSpec2,sizeNodes);Nodes=Nodes(2:end,:);Nodes=Nodes';
fclose(fileID);
fileID = fopen(fileEle);
formatSpec1 = '%d %d %d';sizeFirstLine = [3 1];
FirstLine= fscanf(fileID,formatSpec1,sizeFirstLine);nEle=FirstLine(1);
formatSpec2 = '%d %d %d %d %d';sizeEles=[5,nEle];
Elements = fscanf(fileID,formatSpec2,sizeEles);
Elements=Elements(2:end,:);Elements=Elements';
fclose(fileID);
time_HandleFile=toc(t1);
fprintf('\n Time cost: %.2f',time_HandleFile);
%% Computation grid and density
%%%%%%%% Computation grid
% < ---------------------------- Parameter ---------------------------- > %
xmin=-0.5;dx=0.01;xmax=0.5;
ymin=-0.5;dy=0.01;ymax=0.5;
z0=-0.28;
% < ------------------------------ End  ------------------------------- > %
xgv=xmin:dx:xmax;ygv=ymin:dy:ymax;
nx=length(xgv);ny=length(ygv);
[X2d,Y2d]=meshgrid(xgv,ygv);
%%%%%%%% Density
% < ---------------------------- Parameter ---------------------------- > %
lambda0=1260; % Const
lambda1=[1 2 3]*1000; % Linear
lambda2=[-5 4 -3 3 -4 5]*1000; % Quadratic
lambda3=[8 7 6 5 4 -4 -5 -6 -7 -8]*1000; % Cubic
a000=lambda0;
a100=lambda1(1);a010=lambda1(2);a001=lambda1(3);
a002=lambda2(1);a011=lambda2(2);a020=lambda2(3);
a101=lambda2(4);a110=lambda2(5);a200=lambda2(6);
a003=lambda3(1);a012=lambda3(2);a021=lambda3(3);
a030=lambda3(4);a102=lambda3(5);a111=lambda3(6);
a120=lambda3(7);a201=lambda3(8);a210=lambda3(9);a300=lambda3(10);
rho=@(x,y,z) a000+a100*x+a010*y+a001*z+...
    a002*z.^2+a011*y.*z+a020*y.^2+...
    a101*x.*z+a110*x.*y+a200*x.^2+...
    a003*z.^3+a012*y.*z.^2+a021*y.^2.*z+...
    a030*y.^3+a102*x.*z.^2+a111*x.*y.*z+...
    a120*x.*y.^2+a201*z.*x.^2+a210*y.*x.^2+a300*x.^3;
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Model parameters...');
fprintf('\n Bennu surface model: Number of vertices: %d; Number of faces: %d.',nVert,nFace);
fprintf('\n Bennu tetrahedral mesh model: Number of nodes: %d; Number of elements: %d.',nNode,nEle);
fprintf('\n Density: %s kg/m^3',func2str(rho));
fprintf('\n The external cuboid of the model: X:[%.2f,%.2f] km; Y:[%.2f,%.2f] km; Z:[%.2f,%.2f] km',...
    X1,X2,Y1,Y2,Z1,Z2);
fprintf('\n Computation plane altitude: %.2f km',z0);
fprintf('\n dx: %.2f km; dy: %.2f km; nx: %d; ny: %d; total computation points: %d',dx,dy,nx,ny,nx*ny);
%% Start computaion
%%%%%%%% Analytical solution of RenZY_2017_SG GraPly
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Analytical solution of RenZY_2017_SG GraPly...');
%%%% load result
matFileRenZY=sprintf('Bennu_Variable_%d_%d_z0%g_dx%g_dy%g_gT_RenZY',...
    nVert,nFace,z0,dx,dy);
matFileRenZY = replace(matFileRenZY,{'+','-','.'},{'p','m','d'});
load(matFileRenZY);
%%%% Print the total number of computation points and time cost
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_RenZY);
fprintf('\n Number of facets computed per second: %d faces/sec',nx*ny*nFace/time_RenZY);
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
    =phi_xyz_PolyhedronVI_Parker(xgv,ygv,z0,...
    Nodes,Elements,rho,NT,Mx,My,Ntheta,Nk,eps_f,eps_b,R0,R1,typeNUFFT);
time_polyParker=toc(t2);
%%%% Print algorithm parameters and time cost
fprintf('\n NUFFT-Parker algorithm parameter:');
fprintf('\n kappa: %.2f; zmin: %.4f km; zmax: %.4f km; zref: %.4f km; chi: %.2f; N_T: %d',kappa,Z1,Z2,zref,chi,NT);
fprintf('\n (Mx=My): %d; Ntheta: %d; Nk: %d;eps: %g; R0: %d; R1: %d',...
    Mx,Ntheta,Nk,eps_f,R0,R1);
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_polyParker);
fprintf('\n Number of facets computed per second: %d faces/sec',nx*ny*nFace/time_polyParker);
fprintf('\n Number of elements computed per second: %d tetrahedrons/sec',nx*ny*nEle/time_polyParker);
%% Mapping
%%%%%%%% RenZY_2017_SG GraPly as reference, to observe the error of Fourier method
del_gft_gx=gx_gft-gx_RenZY;del_gft_gy=gy_gft-gy_RenZY;del_gft_gz=gz_gft-gz_RenZY;
del_gft_Txx=Txx_gft-Txx_RenZY;del_gft_Tyy=Tyy_gft-Tyy_RenZY;del_gft_Tzz=Tzz_gft-Tzz_RenZY;
del_gft_Txy=Txy_gft-Txy_RenZY;del_gft_Txz=Txz_gft-Txz_RenZY;del_gft_Tyz=Tyz_gft-Tyz_RenZY;
Datas_RenZY_Fourier={[],[],[];
    gx_RenZY,gx_gft,del_gft_gx;
    gy_RenZY,gy_gft,del_gft_gy;
    gz_RenZY,gz_gft,del_gft_gz;
    Txx_RenZY,Txx_gft,del_gft_Txx;
    Tyy_RenZY,Tyy_gft,del_gft_Tyy;
    Tzz_RenZY,Tzz_gft,del_gft_Tzz;
    Txy_RenZY,Txy_gft,del_gft_Txy;
    Txz_RenZY,Txz_gft,del_gft_Txz;
    Tyz_RenZY,Tyz_gft,del_gft_Tyz
    };
% < ---------------------------- Parameter ---------------------------- > %
fs1=10;fs2=12;
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
Datas=Datas_RenZY_Fourier;
pickGrav=[2,3,4];
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
%%%%%%%% Difference between Fourier and RenZY_2017_SG GraPly solutions
stas_gft_gx=sta_min_max_mean_rms(gx_gft,gx_RenZY,cutEPS_dVg);
stas_gft_gy=sta_min_max_mean_rms(gy_gft,gy_RenZY,cutEPS_dVg);
stas_gft_gz=sta_min_max_mean_rms(gz_gft,gz_RenZY,cutEPS_dVg);
stas_gft_Txx=sta_min_max_mean_rms(Txx_gft,Txx_RenZY,cutEPS_T);
stas_gft_Tyy=sta_min_max_mean_rms(Tyy_gft,Tyy_RenZY,cutEPS_T);
stas_gft_Tzz=sta_min_max_mean_rms(Tzz_gft,Tzz_RenZY,cutEPS_T);
stas_gft_Txy=sta_min_max_mean_rms(Txy_gft,Txy_RenZY,cutEPS_T);
stas_gft_Txz=sta_min_max_mean_rms(Txz_gft,Txz_RenZY,cutEPS_T);
stas_gft_Tyz=sta_min_max_mean_rms(Tyz_gft,Tyz_RenZY,cutEPS_T);
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Difference between Fourier and RenZY_2017_SG GraPly solutions: ');
fprintf('\n min(ref) & max(ref) & mean(ref) & rms(ref) & min(err) & max(err) & mean(err) & rms(err) & E2')
fprintf('\n $g_x$ ');fprintf('& %.2e ',stas_gft_gx(pickErrInfor));fprintf(' \\\\');
fprintf('\n $g_y$ ');fprintf('& %.2e ',stas_gft_gy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $g_z$ ');fprintf('& %.2e ',stas_gft_gz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xx}$ ');fprintf('& %.2e ',stas_gft_Txx(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{yy}$ ');fprintf('& %.2e ',stas_gft_Tyy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{zz}$ ');fprintf('& %.2e ',stas_gft_Tzz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xy}$ ');fprintf('& %.2e ',stas_gft_Txy(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{xz}$ ');fprintf('& %.2e ',stas_gft_Txz(pickErrInfor));fprintf(' \\\\');
fprintf('\n $T_{yz}$ ');fprintf('& %.2e ',stas_gft_Tyz(pickErrInfor));fprintf(' \\\\');















