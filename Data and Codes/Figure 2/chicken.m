clc
clear
%

load ../Data/OTUTables_IC.mat

xs_b = cecumOTU;
[nSamp nOb] = size(xs_b);
nCom = min(nSamp,nOb)
filen = strcat('ModelFits/chicken_',num2str(nCom),'.mat')
% 
ZC     = -0.1*rand(nSamp,nCom);
thetB  = 0.1*rand(nCom,nOb);
load(filen)
% hyperparameters
etaZ = 0.0025;etaT = 0.0025;
grdnorm = 1;iter = 1;ctf_grad = 0.01;
while grdnorm > ctf_grad
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
    deltB = xs_b-QB;

    % gradients
    grthetB = [ZC]'*deltB;
    grz     = deltB*thetB';

    % update the variables
    ZC     = ZC - etaZ*grz;
    thetB  = thetB - etaT*grthetB;
    % errors
    grdnorm = norm(grz)/norm(ZC) + norm(grthetB)/norm(thetB);

    % Output
    if mod(iter,500) == 0
        save(filen,'xs_b','ZC','thetB')
        grdnorm
    end
    iter = iter + 1;
end
plot(xs_b,QB,'o')



