clc
clear
%

clear
K = [1 2 4 8 16 32 64 157];
for k=1:length(K)
    flx = strcat('ModelFits/cows_',num2str(K(k)));
    load(flx)
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);

    JS  = 0.5*(xs_b + QB);
    trm = xs_b.*log(xs_b./JS) + QB.*log(QB./JS);
    trm = nansum(trm');
    mJenShan(k) = mean(trm);
    sJenShan(k) = std(trm);

    trm = bc_pairs(xs_b,QB);
    mBC(k) = mean(trm);
    sBC(k) = std(trm);

end
subplot(1,3,1)
yyaxis left
plot(K,mJenShan,'k')
ylim([0 0.2])
hold on
scatter(K,mJenShan,50,'filled','k')
set(gca,'XScale','log')
yyaxis right
plot(K,mBC,'r')
scatter(K,mBC,50,'filled','r')
(mJenShan-min(mJenShan))/(max(mJenShan)-min(mJenShan))
(mBC-min(mBC))/(max(mBC)-min(mBC))


clear
K = [1 2 4 8 16 32 64 156];
for k=1:length(K)
    flx = strcat('ModelFits/chicken_',num2str(K(k)));
    load(flx)
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);

    JS  = 0.5*(xs_b + QB);
    trm = xs_b.*log(xs_b./JS) + QB.*log(QB./JS);
    trm = nansum(trm');
    mJenShan(k) = mean(trm);
    sJenShan(k) = std(trm);

    trm = bc_pairs(xs_b,QB);
    mBC(k) = mean(trm);
    sBC(k) = std(trm);
end
subplot(1,3,2)
yyaxis left
plot(K,mJenShan,'k')
ylim([0 0.3])
hold on
scatter(K,mJenShan,50,'filled','k')
set(gca,'XScale','log')
yyaxis right
plot(K,mBC,'r')
scatter(K,mBC,50,'filled','r')
ylim([0 0.5])

(mJenShan-min(mJenShan))/(max(mJenShan)-min(mJenShan))
(mBC-min(mBC))/(max(mBC)-min(mBC))

clear
K = [1 2 4 8 16 32 64 95];
for k=1:length(K)
    flx = strcat('ModelFits/human_',num2str(K(k)));
    load(flx)
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);

    JS  = 0.5*(xs_b + QB);
    trm = xs_b.*log(xs_b./JS) + QB.*log(QB./JS);
    trm = nansum(trm');
    mJenShan(k) = mean(trm);
    sJenShan(k) = std(trm);

    trm = bc_pairs(xs_b,QB);
    mBC(k) = mean(trm);
    sBC(k) = std(trm);

end
subplot(1,3,3)
yyaxis left
plot(K,mJenShan,'k')
ylim([0 0.3])
hold on
scatter(K,mJenShan,50,'filled','k')
yyaxis right
plot(K,mBC,'r')
hold on

scatter(K,mBC,50,'filled','r')
set(gca,'XScale','log')
ylim([0 0.5])

(mJenShan-min(mJenShan))/(max(mJenShan)-min(mJenShan))
(mBC-min(mBC))/(max(mBC)-min(mBC))