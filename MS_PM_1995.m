%% P.M. 1995 algorithm

close all
clear all
path = 'D:\MicroStates\MS_data\';
filename = 'BAP-Brisova-avg-19ele-GFPpeaks.txt';

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
%% Explained variance

%variance of data
var_data = trace(M*M')/(Nt*(Ns-1));
expl = (1 - conv/var_data)*100;
%% Visualize
figure(1)
suptitle('P.M. Microstates') %visualize topomaps
for i = 1:size(MS,2)
    subplot(2,2,i)
    topoplot(MS(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end

figure(2)
plot(expl)
title('Explained variance over iteration')
xlabel('Iteration')
ylabel('Eplained variance [%]')

disp(['explained variance by PM model: ',num2str(expl(end))]);
    