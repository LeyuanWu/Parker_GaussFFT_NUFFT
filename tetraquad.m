function [Xs,Ys,Zs,Ws]=tetraquad(Ras,Rbs,Rcs,Rds,N)
% -------------------------------------------------------------------------
% M-file to calculate the N^3 Gaussian nodes and weights for a bunch of
% tetrahedrons simultaneously
% -------------------------------------------------------------------------
% ---- input ---- %
% Ras,Rbs,Rcs,Rds: Coordinates of the four vectices for each tetrahedron (each row)
% N: order of the Gaussian quadrature
% ---- output ---- %
% Xs,Ys,Zs,Ws: the N^3 Gaussian nodes and weights
% -------------------------------------------------------------------------
% ********************
% Modified from
% ********************
% tetraquad.m - Gaussian Quadrature for a tetrahedron
% Reference:
% Greg von Winckel (2021). Gauss Quadrature for Tetrahedra 
% (https://www.mathworks.com/matlabcentral/fileexchange/9389-gauss-quadrature-for-tetrahedra)
% MATLAB Central File Exchange. Retrieved June 5, 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[q1,w1]=rquad(N,2); [q2,w2]=rquad(N,1);  [q3,w3]=rquad(N,0);
[q1,q2,q3]=meshgrid(q1,q2,q3);
q1=q1(:);   q2=q2(:);       q3=q3(:);
x=1-q1;     y=(1-q2).*q1;   z=q1.*q2.*q3;
w=reshape(reshape(w2*w1',N^2,1)*w3',N^3,1);

% c=[1 0 0 0;-1 1 0 0;-1 0 1 0;-1 0 0 1]*Rs;
% Ws=abs(det(c(2:4,:)))*w;

x=transpose(x);y=transpose(y);z=transpose(z);w=transpose(w);
temp1=(Rbs(:,1)-Ras(:,1)).*...
    ((Rcs(:,2)-Ras(:,2)).*(Rds(:,3)-Ras(:,3))-(Rds(:,2)-Ras(:,2)).*(Rcs(:,3)-Ras(:,3)));
temp2=-(Rbs(:,2)-Ras(:,2)).*...
    ((Rcs(:,1)-Ras(:,1)).*(Rds(:,3)-Ras(:,3))-(Rds(:,1)-Ras(:,1)).*(Rcs(:,3)-Ras(:,3)));
temp3=(Rbs(:,3)-Ras(:,3)).*...
    ((Rcs(:,1)-Ras(:,1)).*(Rds(:,2)-Ras(:,2))-(Rds(:,1)-Ras(:,1)).*(Rcs(:,2)-Ras(:,2)));
Ws=abs(temp1+temp2+temp3)*w;

% Change of coordinates 
% XYZ=[ones(N^3,1) x y z]*c; Xs=XYZ(:,1); Ys=XYZ(:,2); Zs=XYZ(:,3);

Xs=repmat(Ras(:,1),1,N^3)+(Rbs(:,1)-Ras(:,1))*x+(Rcs(:,1)-Ras(:,1))*y+(Rds(:,1)-Ras(:,1))*z;
Ys=repmat(Ras(:,2),1,N^3)+(Rbs(:,2)-Ras(:,2))*x+(Rcs(:,2)-Ras(:,2))*y+(Rds(:,2)-Ras(:,2))*z;
Zs=repmat(Ras(:,3),1,N^3)+(Rbs(:,3)-Ras(:,3))*x+(Rcs(:,3)-Ras(:,3))*y+(Rds(:,3)-Ras(:,3))*z;
 
function [x,w]=rquad(N,k)

k1=k+1; k2=k+2; 
n=1:N;  nnk=2*n+k;
A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n);
B1=4*k1/(k2*k2*(k+3)); nk=n+k; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X));    

% Grid points
x=(X+1)/2; 

% Quadrature weights
w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;