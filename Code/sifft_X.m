function F_X=sifft_X(F,xi,alpha)
% -------------------------------------------------------------------------
% M-file to calculate the 1D inverse shift FFT along each row of a matrix F
% -------------------------------------------------------------------------
% ---- input ---- %
% xi: frequency domain shift parameter
% alpha: space domain shift parameter 
% ---- output ---- %
% F_X: 1D inverse shift Fast Fourier Transform of F along each row 
% -------------------------------------------------------------------------
% Reference:
% Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of 
% potential fields: Geophysics, 79, no. 5, G59-G68
% -------------------------------------------------------------------------
F=transpose(F);
F_Y=sifft_Y(F,xi,alpha);
F_X=transpose(F_Y);
