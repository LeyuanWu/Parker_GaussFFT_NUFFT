function stas = sta_min_max_mean_rms(X_cal,X_ref,cutEPS)
% -------------------------------------------------------------------------
% M-file to calculate statistics of a numerical algorithm compared to the
% referenced algorithm
% -------------------------------------------------------------------------
% ---- input ---- %
% X_cal: The numerical algorithm
% X_ref: The referenced algorithm
% cutEPS: % if X_cal < cutEPS*max(X_ref), set eps as nan
% ---- output ---- %
% stas: [min(X_ref),max(X_ref),mean(X_ref),rms(X_ref),
%        min(dif_X),max(dif_X),mean(dif_X),rms(dif_X),
%        min(eps_X),max(eps_X),mean(eps_X),rms(eps_X),
%        E2]
% -------------------------------------------------------------------------
%%%%%%%% Difference and Relative error
dif_X=X_cal-X_ref;
eps_X=(X_cal-X_ref)./X_ref*100;
eps_X(abs(X_ref/max(abs(X_ref(:))))<cutEPS)=nan;
X_ref=X_ref(~isnan(X_ref));
dif_X=dif_X(~isnan(dif_X));
eps_X=eps_X(~isnan(eps_X));
%%%%%%%% Statistics
stas=zeros(1,13); % Row vector easy for table
stas(1)=min(X_ref);
stas(2)=max(X_ref);
stas(3)=mean(X_ref);
stas(4)=sqrt(mean(X_ref.*X_ref));
stas(5)=min(dif_X);
stas(6)=max(dif_X);
stas(7)=mean(dif_X);
stas(8)=sqrt(mean(dif_X.*dif_X));
stas(9)=min(eps_X);
stas(10)=max(eps_X);
stas(11)=mean(eps_X);
stas(12)=sqrt(mean(eps_X.*eps_X));
stas(13)=stas(8)/stas(4);% E2 error
end

