% prove that there is no ENSO in the NoENSO run
addpath(genpath('/gpfs_backup/slarson_data/ktmcmoni/matlab_functions/'))

% load data
ENSO_B=ncread('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/nino34.jan0400-dec0721.anom.nc','SST');
ENSO_MDeq=ncread('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/nino34.jan0685-dec0990.anom.nc','SST');

IOD_B=ncread('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_B_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0400-dec0721.anom.nc','SST');
IOD_MDeq=ncread('/gpfs_backup/slarson_data/slarson/CCSM4/ccsm4_MDeqPac_0.9x1.25_gx1v6/postP/SST/slarson/iod.adj.jan0685-dec0990.anom.nc','SST');

ENSO_B=squeeze(ENSO_B(1,1,end-3671:end));
ENSO_MDeq=squeeze(ENSO_MDeq);

IOD_B=squeeze(IOD_B(1,1,end-3671:end));
IOD_MDeq=squeeze(IOD_MDeq);

[PSI,LAMBDA]=sleptap(length(ENSO_B));
[F1,S1]=mspec((1/12),ENSO_B,PSI,'detrend');
[F2,S2]=mspec((1/12),ENSO_MDeq,PSI,'detrend');
[F3,S3]=mspec((1/12),IOD_B,PSI,'detrend');
[F4,S4]=mspec((1/12),IOD_MDeq,PSI,'detrend');

% plot
figure
set(gcf, 'Position',  [100, 100, 1000, 500])
subplot(1,2,1)
hold on
p1=plot(.5*F1/pi,S1,'k','linewidth',2)
p2=plot(.5*F2/pi,S2,'b','linewidth',2)
uistack(p2,'bottom')
xlabel('Period (years)')
ylabel('PSD (^oC^2 cpy^{-1})')
xlim([0.04 1])
legend([p1 p2],'CTRL','NoENSO','Location','northeast')
xticks([0.04 .1 0.2 0.5 1])
xticklabels({'25','10','5','2','1'})
set(gca,'XScale','log','fontsize',24)
grid on
set(gca,'fontsize',24)
text(0.041,2.8,'a)','fontsize',24)

subplot(1,2,2)
hold on
p1=plot(.5*F3/pi,S3,'k','linewidth',2)
p2=plot(.5*F4/pi,S4,'b','linewidth',2)
uistack(p2,'bottom')
xlabel('Period (years)')
ylabel('PSD (^oC^2 cpy^{-1})')
xlim([0.04 1])
legend([p1 p2],'CTRL','NoENSO','Location','northeast')
xticks([0.04 .1 0.2 0.5 1])
xticklabels({'25','10','5','2','1'})
set(gca,'XScale','log','fontsize',24)
grid on
text(0.041,0.67,'b)','fontsize',24)
set(gca,'fontsize',24)

print -depsc ../analysis/figures/ENSO_spectra.eps

