clear;
fprintf('\n ###################################################### ');
fprintf('\n %s',datetime);
%% Read SRTMv4.1 topography data
%%%%%%%% Himalaya region
lon1=86;lon2=88;lat1=27;lat2=29;
res_deg=1/1200; % Data resolution
%%%%%%%% .tif file name
indLon = floor((lon1+180.0)/5.0)+1;
indLat = floor((60-lat2)/5.0)+1;
srtmFileName=sprintf('srtm_%02d_%02d.tif',indLon,indLat);
%%%%%%%% Read SRTMv4.1 .tif topography data
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Read SRTMv4.1 .tif topography data: %s...',srtmFileName);
t1=tic;
I = GEOTIFF_READ(srtmFileName);
time_srtmRead=toc(t1);
fprintf('\n Complete reading. Time cost: %.2f',time_srtmRead);
%%%%%%%% Grid registered
Lon_deg5=I.x;Lon_deg5=Lon_deg5+res_deg/2;Lon_deg5=Lon_deg5';
Lat_deg5=I.y;Lat_deg5=Lat_deg5-res_deg/2;Lat_deg5=flipud(Lat_deg5');
nLon_deg5=length(Lon_deg5);nLat_deg5=length(Lat_deg5);
[LON_deg5,LAT_deg5]=meshgrid(Lon_deg5,Lat_deg5);
HDEM_m_deg5=I.z;HDEM_m_deg5=flipud(HDEM_m_deg5);
HDEM_m_deg5=double(HDEM_m_deg5); % Convert short to double
%%%%%%%% Extract the part of interest
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Extract the part of interest: %.2f < lon < %.2f, %.2f < lat <%.2f',lon1,lon2,lat1,lat2);
ind_lon1=find(Lon_deg5>=lon1,1,'first');ind_lon2=find(Lon_deg5<=lon2,1,'last');
ind_lat1=find(Lat_deg5>=lat1,1,'first');ind_lat2=find(Lat_deg5<=lat2,1,'last');
Lon=Lon_deg5(ind_lon1:ind_lon2);Lat=Lat_deg5(ind_lat1:ind_lat2);
lonmin=min(Lon);lonmax=max(Lon);dlon=Lon(2)-Lon(1);
latmin=min(Lat);latmax=max(Lat);dlat=Lat(2)-Lat(1);
nLon=length(Lon);nLat=length(Lat);
[LON,LAT]=meshgrid(Lon,Lat);
HDEM_m=HDEM_m_deg5(ind_lat1:ind_lat2,ind_lon1:ind_lon2);
HDEM_km=HDEM_m/1e3;
fprintf('\n Finished extraction, grid size: %d * %d',nLon,nLat);
fprintf(sprintf('\n Data resolution (geographic coordinates): dlon=%g deg, dlat=%g deg',dlon,dlat));
%% Transform to Cartesian coordinates under plane approximation
% < ---------------------------- Parameter ---------------------------- > %
R=6378.1370; % WGS84 Earth Radius in km
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Transform to Cartesian coordinates under plane approximation. ');
fprintf('\n Earth radius: R=%.4f km',R);
%%%%%%%% Center latitude and longitude as the origin (x,y)=(0,0)
dEst_deg=R*cosd(mean(Lat))*pi/180;
dNor_deg=R*pi/180;
Est=(Lon-mean(Lon))*dEst_deg;Estmin=min(Est);Estmax=max(Est);dEst=Est(2)-Est(1);
Nor=(Lat-mean(Lat))*dNor_deg;Normin=min(Nor);Normax=max(Nor);dNor=Nor(2)-Nor(1);
nEst=length(Est);nNor=length(Nor);
[EST,NOR]=meshgrid(Est,Nor);
fprintf(sprintf('\n Cartesian coordinate DEM grid size: %d * %d ',nEst,nNor));
fprintf(sprintf('\n Data resolution (Cartesian coordinates): dEst=%g km, dNor=%g km',dEst,dNor));
%% Model parameters
%%%%%%%% Find the profile passing the peak of the DEM for comparison
[~,col] = find(HDEM_km==max(HDEM_km(:)));
Hmax=max(HDEM_km(:));Hmin=min(HDEM_km(:));
% < ---------------------------- Parameter ---------------------------- > %
rho0=2670;% Density in kg/m^3
href=0; % Reference height
hCal=ceil((Hmax+3*dNor)/0.1)*0.1; % Calculation height. 3 times the grid size above the peak and a multiple of 100 meters
eswn='ne'; % Triangulation direction. 
%            'ne': southwest-northeast triangulation
%            'nw': southeast-northwest triangulation
% < ------------------------------ End  ------------------------------- > %
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Model parameters...');
fprintf('\n Density: %.2f kg/m^3.',rho0);
fprintf('\n Reference height: %.4f km.',href);
fprintf('\n Topographic lowest point: %.4f km; highest point: %.4f km; Calculation height: %.4f km',Hmin,Hmax,hCal);
fprintf('\n Triangulation direction: %s.',eswn);
%% Start computaion
%%%%%%%% Werner_1997_CMDA
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Analytical solution of Werner_1997_CMDA...');
%%%% load result
matFileWerner=sprintf('srtm_%02d_%02d_lon_%+.1f_%+.1f_lat_%+.1f_%+.1f_%s_Werner'...
    ,indLon,indLat,lon1,lon2,lat1,lat2,eswn);
matFileWerner = replace(matFileWerner,{'+','-','.'},{'p','m','d'});
load(matFileWerner);
%%%% Print the total number of computation points and time cost
nP=length(Lat); % Computation carried out only on a profile for space-domain solutions
fprintf('\n Number of computation points: %d',nP);
fprintf('\n Time cost: %.2f sec',time_dVgT_Werner);
%%%%%%%% Prism summation
fprintf('\n ------------------------------------------------------ ');
fprintf('\n Prism summation...');
%%%% load result
matFilePS=sprintf('srtm_%02d_%02d_lon_%+.1f_%+.1f_lat_%+.1f_%+.1f_%s_PrismSum'...
    ,indLon,indLat,lon1,lon2,lat1,lat2,eswn);
matFilePS = replace(matFilePS,{'+','-','.'},{'p','m','d'});
    load(matFilePS);
%%%% Print the total number of computation points and time cost
fprintf('\n Number of computation points: %d',nP);
fprintf('\n Time cost: %.2f sec',time_PS);
%%%%%%%% NUFFT-based Fourier-Parker method
%%%% set x as north, y as east, z downward
xgv=Nor;ygv=Est;dx=dNor;dy=dEst;
nx=length(xgv);ny=length(ygv);
[X2d,Y2d]=meshgrid(xgv,ygv);
z0=-hCal;
%%%% Parameters related to the convergence behavior of Parker's formula
Z1=-max(HDEM_km(:));Z2=-min(HDEM_km(:));
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
    =phi_xyz_DEMSI_Parker(xgv,ygv,z0,...
    HDEM_km,href,nEst,nNor,Estmin,Estmax,Normin,Normax,eswn,rho0,...
    NT,Mx,My,Ntheta,Nk,eps_f,eps_b,R0,R1,typeNUFFT);
time_polyParker=toc(t2);
%%%% Print algorithm parameters and time cost
fprintf('\n NUFFT-Parker algorithm parameter:');
fprintf('\n kappa: %.2f; zmin: %.4f km; zmax: %.4f km; zref: %.4f km; chi: %.2f; N_T: %d',kappa,Z1,Z2,zref,chi,NT);
fprintf('\n (Mx=My): %d; Ntheta: %d; Nk: %d;eps: %g; R0: %d; R1: %d',...
    Mx,Ntheta,Nk,eps_f,R0,R1);
fprintf('\n Number of computation points: %d',nx*ny);
fprintf('\n Time cost: %.2f sec',time_polyParker);
%% Mapping
% < ---------------------------- Parameter ---------------------------- > %
fs=10;MS=[7,4,2];
% < ------------------------------ End  ------------------------------- > %
%%%%%%%% Results on the profile
Xs=Lat;
%%%% Werner_1997_CMDA
Werner_POTVEC={dV,gx,gy,gz};
Werner_TEN={Txx,Tyy,Tzz,Txy,Txz,Tyz};
%%%% Prism summation
PS_POTVEC={dV_PS,gx_PS,gy_PS,gz_PS};
PS_TEN={Txx_PS,Tyy_PS,Tzz_PS,Txy_PS,Txz_PS,Tyz_PS};
%%%% NUFFT Parker results, notice the transpose between (Est,Nor) and (X,Y)
dV_gft_line=dV_gft(col,:)';
gx_gft_line=gx_gft(col,:)';gy_gft_line=gy_gft(col,:)';gz_gft_line=gz_gft(col,:)';
Txx_gft_line=Txx_gft(col,:)';Tyy_gft_line=Tyy_gft(col,:)';Tzz_gft_line=Tzz_gft(col,:)';
Txy_gft_line=Txy_gft(col,:)';Txz_gft_line=Txz_gft(col,:)';Tyz_gft_line=Tyz_gft(col,:)';
Parker_POTVEC={dV_gft_line,gx_gft_line,gy_gft_line,gz_gft_line};
Parker_TEN={Txx_gft_line,Tyy_gft_line,Tzz_gft_line,Txy_gft_line,Txz_gft_line,Tyz_gft_line};
%%%%%%%% Potential and vector gravity fields
figure1= figure('Color',[1 1 1]);
YLabelS={'$V$ (m$^2$/s$^2$)','$\delta(V)$ (m$^2$/s$^2$)';
    '$g_x$ (mGal)','$\delta(g_x)$ (mGal)';
    '$g_y$ (mGal)','$\delta(g_y)$ (mGal)';
    '$g_z$ (mGal)','$\delta(g_z)$ (mGal)';
    };
%%%% Comparison
for i=1:1:4
    Y1=Werner_POTVEC{i};Y2=PS_POTVEC{i};Y3=Parker_POTVEC{i};
    subplot1=subplot(4,2,2*i-1,'Parent',figure1);
    Ys=horzcat(Y1,Y2,Y3);
    dispNames={'Polyhedron','PS','Fourier'};
    myplot(subplot1,Xs,Ys,fs,MS,dispNames,0);
    xlabel('Latitude');ylabel(YLabelS{i,1},'Interpreter','latex');
    % Werner_1997_CMDA as reference
    subplot1=subplot(4,2,2*i,'Parent',figure1);
    Ys=horzcat(Y2-Y1,Y3-Y1);
    dispNames={'PS','Fourier'};
    myplot(subplot1,Xs,Ys,fs,MS,dispNames,1);
    xlabel('Latitude');ylabel(YLabelS{i,2},'Interpreter','latex');
end
%%%%%%%% Gravity Gradient Tensor
figure2= figure('Color',[1 1 1]);
YLabelS={'$T_{xx}$ (Eotvos)','$\delta(T_{xx})$ (Eotvos)';
    '$T_{yy}$ (Eotvos)','$\delta(T_{yy})$ (Eotvos)';
    '$T_{zz}$ (Eotvos)','$\delta(T_{zz})$ (Eotvos)';
    '$T_{xy}$ (Eotvos)','$\delta(T_{xy})$ (Eotvos)';
    '$T_{xz}$ (Eotvos)','$\delta(T_{xz})$ (Eotvos)';
    '$T_{yz}$ (Eotvos)','$\delta(T_{yz})$ (Eotvos)'};
%%%% Comparison
for i=1:1:6
    Y1=Werner_TEN{i};Y2=PS_TEN{i};Y3=Parker_TEN{i};
    subplot1=subplot(6,2,2*i-1,'Parent',figure2);
    Ys=horzcat(Y1,Y2,Y3);
    dispNames={'Polyhedron','PS','Fourier'};
    myplot(subplot1,Xs,Ys,fs,MS,dispNames,0);
    xlabel('Latitude');ylabel(YLabelS{i,1},'Interpreter','latex');
    % Werner_1997_CMDA as reference
    subplot1=subplot(6,2,2*i,'Parent',figure2);
    Ys=horzcat(Y2-Y1,Y3-Y1);
    dispNames={'PS','Fourier'};
    myplot(subplot1,Xs,Ys,fs,MS,dispNames,1);
    xlabel('Latitude');ylabel(YLabelS{i,2},'Interpreter','latex');
end
%% Print information on screen (also for table)
% < ---------------------------- Parameter ---------------------------- > %
cutEPS_dVg=1e-3; % if value < cutEPS_dVg*max, set eps as nan
cutEPS_T=1e-3; % if value < cutEPS_T*max, set eps as nan
% < ------------------------------ End  ------------------------------- > %
% < ---------------------------- Parameter ---------------------------- > %
pickErrInfor=[1:8,13];
% < ------------------------------ End  ------------------------------- > %
%%%%%%%% Difference between Fourier and Werner_1997_CMDA solutions
stas_gft_dV=sta_min_max_mean_rms(dV_gft_line,dV,cutEPS_dVg);
stas_gft_gx=sta_min_max_mean_rms(gx_gft_line,gx,cutEPS_dVg);
stas_gft_gy=sta_min_max_mean_rms(gy_gft_line,gy,cutEPS_dVg);
stas_gft_gz=sta_min_max_mean_rms(gz_gft_line,gz,cutEPS_dVg);
stas_gft_Txx=sta_min_max_mean_rms(Txx_gft_line,Txx,cutEPS_T);
stas_gft_Tyy=sta_min_max_mean_rms(Tyy_gft_line,Tyy,cutEPS_T);
stas_gft_Tzz=sta_min_max_mean_rms(Tzz_gft_line,Tzz,cutEPS_T);
stas_gft_Txy=sta_min_max_mean_rms(Txy_gft_line,Txy,cutEPS_T);
stas_gft_Txz=sta_min_max_mean_rms(Txz_gft_line,Txz,cutEPS_T);
stas_gft_Tyz=sta_min_max_mean_rms(Tyz_gft_line,Tyz,cutEPS_T);
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
