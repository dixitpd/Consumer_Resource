clc
clear
%

load ../Data/three_kingdoms_cleaned_up

% 
%% Only abundant OTUs in species 1 of cows
sps1 = find(cow_species==1);
ctf_otu      = 0.001;% based on Brian's paper
bacteria_xs  = bacteria_xs(sps1,:);
mn_bact      = mean(bacteria_xs);
goodb        = find(mn_bact > ctf_otu);
restb = 1 - sum(bacteria_xs(:,goodb)')';
xs_b  = [bacteria_xs(:,goodb) restb];
% 
load ModelFits/test_data
xs_b_train = xs_b;
xs_b_train(test_data,:) = [];
xs_b_test = xs_b(test_data,:);
% 
[nSamp nOb] = size(xs_b_train);
nCom = 16;

filen = strcat('ModelFits/cows_train_test_',num2str(nCom),'.mat')

ZC     = -0.1*rand(nSamp,nCom);
thetB  = 0.1*rand(nCom,nOb);
%load(filen)
% hyperparameters
etaZ = 0.005;etaT = 0.005;
grdnorm = 1;iter = 1;ctf_grad = 0.01;
while grdnorm > ctf_grad
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
    deltB = xs_b_train-QB;

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
        save(filen,'xs_b_train','xs_b_test','ZC','thetB')
        grdnorm
    end
    iter = iter + 1;
end
plot(QB,xs_b_train,'o')
