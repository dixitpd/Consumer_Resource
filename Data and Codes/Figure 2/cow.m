clc
clear
%

load ../Data/three_kingdoms_cleaned_up


%% Only abundant OTUs in species 1 of cows
sps1     = find(cow_species==1);
ctf_otu  = 0.001;% based on Brian's paper
xs_1     = bacteria_xs(sps1,:);
mn_bact  = mean(xs_1);
goodb1   = find(mn_bact > ctf_otu);

sps2     = find(cow_species==2);
ctf_otu  = 0.001;% based on Brian's paper
xs_2     = bacteria_xs(sps2,:);
mn_bact  = mean(xs_2);
goodb2   = find(mn_bact > ctf_otu);
goodb    = union(goodb1,goodb2);
restb = 1 - sum(bacteria_xs(:,goodb)')';
xs_b  = [bacteria_xs(:,goodb) restb];
[nSamp nOb] = size(xs_b);
% nCom = min(nSamp,nOb);
% 
% filen = strcat('ModelFits/cows_',num2str(nCom),'.mat')
% 
% ZC     = -0.1*rand(nSamp,nCom);
% thetB  = 0.1*rand(nCom,nOb);
% %load(filen)
% % hyperparameters
% etaZ = 0.005;etaT = 0.005;
% grdnorm = 1;iter = 1;ctf_grad = 0.01;
% while grdnorm > ctf_grad
%     QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
%     deltB = xs_b-QB;
% 
%     % gradients
%     grthetB = [ZC]'*deltB;
%     grz     = deltB*thetB';
% 
%     % update the variables
%     ZC     = ZC - etaZ*grz;
%     thetB  = thetB - etaT*grthetB;
%     %ZC(ZC > 0) = 0;
% 
%     % errors
%     grdnorm = norm(grz)/norm(ZC) + norm(grthetB)/norm(thetB);
% 
%     % Output
%     if mod(iter,200) == 0
%         save(filen,'xs_b','ZC','thetB')
%         grdnorm
%     end
%     iter = iter + 1;
% end
% plot(QB,xs_b,'o')
% % thetB = thetB - min(thetB')';
% % t = ZC*thetB;
% % [u s v] = svd(t);mx = min(nOb,nCom);
% % s = diag(s);s = s(1:mx);
% % f = s(1:end-1)./s(2:end);
% % loglog(f,'o')
% % mn = mean(f(20:end-20));
% % st = std(f(20:end-20));
% % hold on
% % plot(1:100,(mn+st)*ones(100,1),'r--')
% % plot(1:100,(mn-st)*ones(100,1),'r--')
% % xlim([1 100])
% % 
% % 
% % 
