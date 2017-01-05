%% PCA vs. varimax microstates

close all
clear all
path = 'D:\MicroStates\MS_data\';
filename = 'BAP-Brisova-avg-19ele-GFPpeaks.txt';

M = dlmread([path,filename],' ');

Ns = 19;
M = M(:,1:Ns);
Nt = size(M,1);
H = eye(Ns) - ones(Ns)/Ns;  %linear average reference transformation mat.
M = (M * H)';               %average reference and transpose
M_stand = zscore(M,0,1);    %standardize all topomaps
M_stand = M;

%% Relaxed clusturing approach
M_pca_components = pca(M_stand);       %initial factors from PCA
MS_pca = M_stand * M_pca_components(:,1:4); %project original data into PCA space
[M_varimax_components,T] = rotatefactors(M_pca_components(:,1:4)); %rotate the first four (varimax is default)
MS_varmax = M_stand * M_varimax_components;             %project original data into rotated space

%normalize the microStates
MS_pca = (MS_pca'./repmat(sqrt(diag(MS_pca'*MS_pca)),1,19))';
MS_varmax = (MS_varmax'./repmat(sqrt(diag(MS_varmax'*MS_varmax)),1,19))';


%% PM approach
PM_pca_components = pca(M_stand');
MS_PM_pca = PM_pca_components(:,1:4);   %here the components are the MS
[MS_PM_varimax,T] = rotatefactors(PM_pca_components(:,1:4));

%% Visualize
figure
suptitle('Varimax Microstates - spatial') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_varmax(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('PCA Microstates - spatial') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_pca(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('Varimax Microstates - temporal') %visualize topomaps
for i = 1:size(MS_PM_varimax,2)
    subplot(2,2,i)
    topoplot(MS_PM_varimax(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('PCA Microstates - temporal') %visualize topomaps
for i = 1:size(MS_PM_pca,2)
    subplot(2,2,i)
    topoplot(MS_PM_pca(:,i),'10-20-microst.loc','gridscale',150);
end

%% Mapping the relaxed solution to MS constrains

M_pca_components = M_pca_components(:,1:4);   %pick the first four components

for i = 1:size(M_pca_components,1)
    vec = M_pca_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_pca_labels(i,:) = vec/max(abs(vec));
    
    vec = M_varimax_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_varimax_labels(i,:) = vec/max(abs(vec));
end

%% Explained variance
M_stand = M_stand';
aprox_varimax_norm = M_varimax_labels * MS_varmax';
aprox_pca_norm = M_pca_labels * MS_pca';

%measured reached residual variance
varimax_loadings = diag(M_stand*aprox_varimax_norm');
pca_loadings = diag(M_stand*aprox_pca_norm');
aprox_varimax = aprox_varimax_norm.*repmat(varimax_loadings,1,19);
aprox_pca = aprox_pca_norm.*repmat(pca_loadings,1,19);

res_var_varimax = trace((M_stand - aprox_varimax)*(M_stand - aprox_varimax)')/(Nt*(Ns-1));
res_var_pca = trace((M_stand - aprox_pca)*(M_stand - aprox_pca)')/(Nt*(Ns-1));

%theoretical minimum residual variance
min_res_var_varimax = trace(M_stand*M_stand' - (aprox_varimax_norm*M_stand').^2)/(Nt*(Ns-1));
min_res_var_pca = trace(M_stand*M_stand' - (aprox_pca_norm*M_stand').^2)/(Nt*(Ns-1));

%variance of data
var_data = trace(M_stand*M_stand')/(Nt*(Ns-1));

%theoretical explained variance
theo_explained_varimax = 1 - min_res_var_varimax/var_data;
theo_explained_pca = 1 - min_res_var_pca/var_data;
%numerical explained variance
explained_varimax = 1 - res_var_varimax/var_data;
explained_pca = 1 - res_var_pca/var_data;

disp(['theoretical explained variance by varimax: ',num2str(theo_explained_varimax)]);
disp(['numerical explained variance by varimax: ',num2str(explained_varimax)]);
disp(['theoretical explained variance by pca: ',num2str(theo_explained_pca)]);
disp(['numerical explained variance by pca: ',num2str(explained_pca)]);


%% Additional visualization

M_pca_sq = M_pca_components.^2;   
M_varimax_sq = M_varimax_components.^2;

