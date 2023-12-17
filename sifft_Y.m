function F_Y=sifft_Y(F,eta,beta)
% -------------------------------------------------------------------------
% M-file to calculate the 1D inverse shift FFT along each column of a matrix F
% -------------------------------------------------------------------------
% ---- input ---- %
% eta: frequency domain shift parameter
% beta: space domain shift parameter 
% ---- output ---- %
% F_Y: 1D inverse shift Fast Fourier Transform of F along each column 
% -------------------------------------------------------------------------
% Reference:
% Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of 
% potential fields: Geophysics, 79, no. 5, G59-G68
% -------------------------------------------------------------------------
[ny,nx]=size(F);
ny1=ceil(-ny/2);
ny2=ceil(ny/2-1);
NY=(ny1:1:ny2)';

wy=exp(1i*2*pi/ny);
TR_y=wy.^(beta*(NY+eta));
TR_Y=repmat(TR_y,[1,nx]);
TL_y=wy.^(eta*NY);
TL_Y=repmat(TL_y,[1,nx]);
F1=ifftshift(F.*TR_Y);
X1=ifft(F1);
F_Y=TL_Y.*fftshift(X1);