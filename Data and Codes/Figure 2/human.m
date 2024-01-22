clc
clear
%


load ../Data/Human/combined_met.mat
xs_b = motuncbispecies;
badB = find(mean(xs_b)<0.001);
restb = sum(xs_b(:,badB)');xs_b(:,badB) = [];
xs_b = [xs_b restb'];
[nSamp nOb] = size(xs_b);
nCom = min(nSamp,nOb);
filen = strcat('ModelFits/human_',num2str(nCom),'.mat')

ZC     = -0.1*rand(nSamp,nCom);
thetB  = 0.1*rand(nCom,nOb);
%load(filen)
% hyperparameters
etaZ = 0.005;etaT = 0.005;
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
    %ZC(ZC > 0) = 0;

    % errors
    grdnorm = norm(grz)/norm(ZC) + norm(grthetB)/norm(thetB);

    % Output
    if mod(iter,200) == 0
        save(filen,'xs_b','ZC','thetB')
        grdnorm
    end
    iter = iter + 1;
end

plot(QB,xs_b,'o')


