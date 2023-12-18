function X=sifft2(F,xi,eta,alpha,beta)
% -------------------------------------------------------------------------
% M-file to calculate the 2D inverse shift FFT of a spectrum F (2D matrix)
% -------------------------------------------------------------------------
% ---- input ---- %
% alpha,beta: space domain shift parameter of X along two dimensions
% xi,eta: frequency domain shift parameter of F along two dimensions
% ---- output ---- %
% X: 2D inverse shift Fast Fourier Transform of F 
% -------------------------------------------------------------------------
% Reference:
% Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of 
% potential fields: Geophysics, 79, no. 5, G59-G68
% -------------------------------------------------------------------------
F_Y=sifft_Y(F,eta,beta);
X=sifft_X(F_Y,xi,alpha);
