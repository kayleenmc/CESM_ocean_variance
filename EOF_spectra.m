% we want to take the spectra of the first PC timeseries from each experiment
% run it in matlab because python spectra not great

addpath(genpath('/gpfs_backup/slarson_data/ktmcmoni/matlab_functions/'))
load ../EOFs.mat

if 1 % move units BACK to PC
PC1_B=PC1_B.*std(EOF1_B)
PC1_MD=PC1_MD.*std(EOF1_MD)
PC1_MDeq=PC1_MDeq.*std(EOF1_MDeq)
end

[PSI,LAMBDA]=sleptap(length(PC1_B));
[F1,S1]=mspec((1/12),PC1_B.',PSI,'detrend');
%[F2,S2]=mspec((1/12),PC1_MD.',PSI,'detrend');
[F3,S3]=mspec((1/12),PC1_MDeq.',PSI,'detrend');

[PSI,LAMBDA]=sleptap(length(PC1_MD));
[F2,S2]=mspec((1/12),PC1_MD.',PSI,'detrend');

[ra,rb]=mconf(7,0.95)

if 1
figure
hold on
p1=plot(.5*F1/pi,S1,'k','linewidth',2)
%p2=plot(.5*F2/pi,S2,'r','linewidth',2)
p3=plot(.5*F3/pi,S3,'b','linewidth',2)
%uistack(p2,'bottom')
xlabel('Frequency (cpy)')
ylabel('PSD (PW/cpy)')
legend([p1 p3],'FC','no ENSO','Location','southwest')
set(gca,'XScale','log','YScale','log','fontsize',24)
%title('Power spectra of PC1')
grid on
set(gca,'fontsize',24)
print -dpng ../analysis/figures/EOF_spectra2.png
end

if 1
% make a figure showing the pattern of each EOF, same colors as plot above
figure
hold on
p1=plot(EOF1_B,-35:19,'k','linewidth',2)
p2=plot(-1*EOF1_MDeq,-35:19,'b','linewidth',2)
%p3=plot(EOF1_MD,-35:19,'r','linewidth',2)
xlabel('PW (10^{15} W)')
ylabel('Latitude')
legend([p1 p2],'FC','NoENSO','Location','northwest')
grid on
%title('EOF1 of MHT')
set(gca,'fontsize',24)
print -dpng ../analysis/figures/EOF_IO_single2.png
end

% make a spectra like Sarah's

% make a dummy variable for period
figure
hold on
p1=plot(.5*F1/pi,S1.*(.5*F1/pi),'k','linewidth',2)
%plot(.5*F1/pi,S1*ra,'k:','linewidth',1)
%plot(.5*F1/pi,S1*rb,'k:','linewidth',1)
%p2=plot(.5*F2/pi,S2.*(.5*F2/pi),'r','linewidth',2)
%plot(.5*F2/pi,S2*ra,'r:','linewidth',1)
%plot(.5*F2/pi,S2*rb,'r:','linewidth',1)
p3=plot(.5*F3/pi,S3.*(.5*F3/pi),'b','linewidth',2)
%plot(.5*F3/pi,S3*ra,'b:','linewidth',1)
%plot(.5*F3/pi,S3*rb,'b:','linewidth',1)
%uistack(p2,'bottom')
xlabel('Period (years)')
ylabel('PSD (PW^2/cpy^2)')
xlim([0.02 1])
%ylim([0 0.25])
legend([p1 p3],'FC','no ENSO','Location','northwest')
xticks([0.02 0.04 .1 0.2 0.5 1])
xticklabels({'50','25','10','5','2','1'})
set(gca,'XScale','log','fontsize',24)
%title('PC1 of Indian MHT')
grid on
set(gca,'fontsize',24)
print -dpng ../analysis/figures/EOF_spectra_zoom_preserve_var2.png

figure
hold on
p1=plot(.5*F1/pi,S1,'k','linewidth',2)
%plot(.5*F1/pi,S1*ra,'k:','linewidth',1)
%plot(.5*F1/pi,S1*rb,'k:','linewidth',1)
%p2=plot(.5*F2/pi,S2,'r','linewidth',2)
%plot(.5*F2/pi,S2*ra,'r:','linewidth',1)
%plot(.5*F2/pi,S2*rb,'r:','linewidth',1)
p3=plot(.5*F3/pi,S3,'b','linewidth',2)
%plot(.5*F3/pi,S3*ra,'b:','linewidth',1)
%plot(.5*F3/pi,S3*rb,'b:','linewidth',1)
uistack(p2,'bottom')
xlabel('Period (years)')
ylabel('PSD (PW/cpy^2)')
xlim([0.02 1])
%ylim([0 0.25])
legend([p1 p3],'FC','no ENSO','Location','northwest')
xticks([0.02 0.04 .1 0.2 0.5 1])
xticklabels({'50','25','10','5','2','1'})
set(gca,'XScale','log','fontsize',24)
%title('PC1 of Indian MHT')
grid on
set(gca,'fontsize',24)
print -dpng ../analysis/figures/EOF_spectra_zoom2.png

% subplot of EOF, spectra
figure
set(gcf, 'Position',  [100, 100, 500, 1000])
subplot(2,1,1)
hold on
p1=plot(EOF1_B,-35:19,'k','linewidth',2)
p2=plot(-1*EOF1_MDeq,-35:19,'b','linewidth',2)
%p3=plot(EOF1_MD,-35:19,'r','linewidth',2)
xlabel('PW (10^{15} W)')
ylabel('Latitude')
legend([p1 p2],'CTRL','NoENSO','Location','northwest')
grid on
%title('EOF1 of MHT')
set(gca,'fontsize',24)
%text(-0.39,-38,'a)','fontsize',24)
annotation('textbox',[0.24,0.96,0,0],'string','a)','fontsize',24)

subplot(2,1,2)
hold on
p1=plot(.5*F1/pi,S1,'k','linewidth',2)
p3=plot(.5*F3/pi,S3,'b','linewidth',2)
uistack(p2,'bottom')
xlabel('Period (years)')
ylabel('PSD (PW^2 cpy^{-1})')
xlim([0.0067 1])
legend([p1 p3],'CTRL','NoENSO','Location','northwest')
xticks([0.0067 0.02 0.04 .1 0.2 0.5 1])
xticklabels({'150','50','25','10','5','2','1'})
set(gca,'XScale','log','fontsize',24)
grid on
set(gca,'fontsize',24)
annotation('textbox',[0.24,0.48,0,0],'string','b)','fontsize',24)
%text(0.008,0.002,'b)','fontsize',24)
print -depsc ../analysis/figures/EOF_spectra_zoom_subplots.eps
