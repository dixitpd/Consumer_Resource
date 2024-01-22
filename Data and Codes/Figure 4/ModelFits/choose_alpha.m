clc
clear

alx = [0.8 0.9 0.95 0.975 0.995];
for iter=1:length(alx)
    filen = strcat('cows_metadata_24_',num2str(alx(iter)))
    load(filen)
    kl(iter) = mean(-nansum( (xs_b_train.*log(QB./xs_b_train))' ));
    l2(iter) = norm(ZC*C - md_train);
end

subplot(1,2,1)
hold on
yyaxis left
plot(alx,kl,'k')
scatter(alx,kl,20,'k','filled')
yyaxis right
plot(alx,l2,'r')
scatter(alx,l2,20,'r','filled')


alx = [0.5 0.8 0.9 0.95 0.975 0.995];

for iter=1:length(alx)
    filen = strcat('human_metadata_24_',num2str(alx(iter)))
    load(filen)
    QB = exp(-ZC*thetB);QB = normalize(QB,2,'norm',1);
    QM = exp(-ZC*thetM);
    kl(iter) = mean(-nansum( (xs_b_train.*log(QB./xs_b_train))' ));
    l2(iter) = sum(nansum(metconcs_train.*log(metconcs_train./QM),2)) + sum(sum(QM)) - sum(sum(metconcs_train));
end

subplot(1,2,2)
hold on
yyaxis left
plot(alx,kl,'k')
scatter(alx,kl,20,'k','filled')
yyaxis right
plot(alx,l2,'r')
scatter(alx,l2,20,'r','filled')
