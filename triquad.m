function [Xs,Ys,Zs,Ws,nis]=triquad(Ras,Rbs,Rcs,N)
% -------------------------------------------------------------------------
% M-file to calculate the N^2 Gaussian nodes and weights for a bunch of
% 3D triangles simultaneously
% -------------------------------------------------------------------------
% ---- input ---- %
% Ras,Rbs,Rcs: Coordinates of the three vectices for each triangle (each row)
% N: order of the Gaussian quadrature
% ---- output ---- %
% Xs,Ys,Zs,Ws: the N^2 Gaussian nodes and weights
% nis: the outer normals for each triangle
% -------------------------------------------------------------------------
% ********************
% Modified from
% ********************
% triquad.m - Gaussian Quadrature for a triangular domain
% Reference:
% Greg von Winckel (2021). Gaussian Quadrature for Triangles 
% (https://www.mathworks.com/matlabcentral/fileexchange/9230-gaussian-quadrature-for-triangles)
% MATLAB Central File Exchange. Retrieved June 5, 2021.
nF=size(Ras,1);
Mabs=Rbs-Ras;
Macs=Rcs-Ras;
Us=Mabs./sqrt(sum(Mabs.^2,2));
nTemps=horzcat(Mabs(:,2).*Macs(:,3)-Mabs(:,3).*Macs(:,2),...
    Mabs(:,3).*Macs(:,1)-Mabs(:,1).*Macs(:,3),...
    Mabs(:,1).*Macs(:,2)-Mabs(:,2).*Macs(:,1));
nis=nTemps./sqrt(sum(nTemps.^2,2));
Vs=horzcat(nis(:,2).*Us(:,3)-nis(:,3).*Us(:,2),...
    nis(:,3).*Us(:,1)-nis(:,1).*Us(:,3),...
    nis(:,1).*Us(:,2)-nis(:,2).*Us(:,1));
Ras2d=zeros(nF,2);
Rbs2d=horzcat(sqrt(sum(Mabs.^2,2)),zeros(nF,1));
Rcs2d=horzcat(sum(Macs.*Us,2),sum(Macs.*Vs,2));
[tXs,tYs,Ws]=triquad2d_v2(Ras2d,Rbs2d,Rcs2d,N);
Xs=repmat(Ras(:,1),1,N^2)+repmat(Us(:,1),1,N^2).*tXs+repmat(Vs(:,1),1,N^2).*tYs;
Ys=repmat(Ras(:,2),1,N^2)+repmat(Us(:,2),1,N^2).*tXs+repmat(Vs(:,2),1,N^2).*tYs;
Zs=repmat(Ras(:,3),1,N^2)+repmat(Us(:,3),1,N^2).*tXs+repmat(Vs(:,3),1,N^2).*tYs;
end

function [Xs,Ys,Ws]=triquad2d_v2(Ras2d,Rbs2d,Rcs2d,N)

n=1:N;  nnk=2*n+1; A=[1/3 ones(1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
[V,Xs]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[Xs,I]=sort(diag(Xs)); x=(Xs+1)/2; Wx=ab(1,2)*V(1,I)'.^2/4;

N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
L=zeros(N1,N2);  y0=2;  iter=0;
while max(abs(y-y0))>eps
    L(:,1)=1;    L(:,2)=y;
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
end
t=(1+y)/2;
Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
% cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v;
W=Wx*Wy';
Ws=abs(Rcs2d(:,1).*(Rbs2d(:,2)-Ras2d(:,2))...
    +Rbs2d(:,1).*(Ras2d(:,2)-Rcs2d(:,2))...
    +Ras2d(:,1).*(Rcs2d(:,2)-Rbs2d(:,2)))*transpose(W(:));
[tt,xx]=meshgrid(t,x); yy=tt.*xx;
xx=transpose(xx(:));yy=transpose(yy(:));
Xs=Ras2d(:,1)+(Rcs2d(:,1)-Ras2d(:,1))*xx+(Rbs2d(:,1)-Rcs2d(:,1))*yy;
Ys=Ras2d(:,2)+(Rcs2d(:,2)-Ras2d(:,2))*xx+(Rbs2d(:,2)-Rcs2d(:,2))*yy;
end