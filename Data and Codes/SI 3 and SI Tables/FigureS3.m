clc
clear

clc
clear

alx = [0.8 0.9 0.95 0.975 0.995];
for iter=1:length(alx)
    filen = strcat('../Figure4/ModelFits/cows_metadata_24_',num2str(alx(iter)))
    load(filen)
    kl(iter) = mean(-nansum( (xs_b_train.*log(QB./xs_b_train))' ));
    l2(iter) = norm(ZC*C - md_train);
end

hold on
yyaxis left
plot(alx,kl,'k')
scatter(alx,kl,20,'k','filled')
yyaxis right
plot(alx,l2,'r')
scatter(alx,l2,20,'r','filled')


