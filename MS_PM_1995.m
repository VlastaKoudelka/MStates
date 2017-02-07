%% P.M. 1995 algorithm

close all
clear all
% path = 'D:\MicroStates\MS_data\';
path = 'D:\MicroStates\KubaK\pla008\';
% filename = 'BAP-Brisova-avg-19ele-GFPpeaks.txt';
filename = 'PLA_008_EEG_B_r01604-20150223-121157_Edit Channels 2.DAT';

M = dlmread([path,filename],' ');

Ns = 19;
Nu = 4;
Nt = size(M,1);
sig0 = 0;
sigu = 1; %an initial value of residual MS model variance
eps = 1e-6;
M = M(:,1:Ns);
Nt = size(M,1);
H = eye(Ns) - ones(Ns)/Ns;  %linear average reference transformation mat.
M = (M * H);               %average reference and transpose

%Generate initial labels step 2b) in P.M. 1995
labels = (rand(Nt,Nu));
for t = 1:Nt
    vec = labels(t,:);
    vec(abs(vec)<max(abs(vec)))=0;
    labels(t,:)=ceil(vec);
end
i = 1;
while abs(sig0 - sigu) > eps*sigu %step 6)
    sig0 = sigu;
    
    %step 4)
    for k = 1:Nu        %over all MicroStates
        S{k} = zeros(Ns);
        for t = 1:Nt    %over all time samples
            if labels(t,k) == 1     %accumulate only selected samples
                S{k} = S{k} + M(t,:)'*M(t,:);
            end
        end
        [V,D]=eig(S{k});    %derive eigenvectors
        MS(:,k)=V(:,end);   %the last column eigenvect. has the largest
    end 
   
    %step 5)
    MS_aprox = labels*MS';
    sigu = trace(M*M' - (MS_aprox*M').^2)/(Nt*(Ns-1));   
    conv(i) = sigu;
    
    %step 3)
    dist = (M*MS).^2;
    for t = 1:Nt
        vec = dist(t,:);
        vec(vec<max(vec))=0;
        vec(vec~=0)=1;
        labels(t,:)=vec;        
    end
    i = i+1;
end
%% Relaxed solution
for i = 1:Nt
    c(i,:) = inv(MS'*MS)*MS'*M(i,:)';   %least square fit
end

aprox_MS_relaxed = c*MS';

%% Explained variance

%variance of data
var_data = trace(M*M')/(Nt*(Ns-1));

%residual variance or relaxed model
res_var_pca_relaxed = trace((M - aprox_MS_relaxed)*(M - aprox_MS_relaxed)')/(Nt*(Ns-1));


%explained by restricted and relaxed models
expl_restrict = (1 - conv/var_data)*100;
expl_relaxed = (1 - res_var_pca_relaxed/var_data)*100;

%% Visualize
figure(1)
suptitle('P.M. Microstates') %visualize topomaps
for i = 1:size(MS,2)
    subplot(2,2,i)
    topoplot(MS(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end

figure(2)
plot(expl_restrict)
title('Explained variance over iteration')
xlabel('Iteration')
ylabel('Eplained variance [%]')

disp(['explained variance by PM model: ',num2str(expl_restrict(end))]);
disp(['explained variance by relaxed PM model: ',num2str(expl_relaxed)]);
    