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

%% Relaxed clusturing approach - spatial
M_pca_components = pca(M);       %initial factors from PCA
MS_pca = M * M_pca_components(:,1:4); %project original data into PCA space
[M_varmax_components,T] = rotatefactors(M_pca_components(:,1:4)); %rotate the first four (varimax is default)
MS_varmax = M * M_varmax_components;             %project original data into rotated space

%normalize the microStates
MS_pca_norm = (MS_pca'./repmat(sqrt(diag(MS_pca'*MS_pca)),1,19))';
MS_varmax_norm = (MS_varmax'./repmat(sqrt(diag(MS_varmax'*MS_varmax)),1,19))';


%% PM approach - temporal
PM_pca_components = pca(M');
MS_PM_pca = PM_pca_components(:,1:4);   %here the components are the MS
[MS_PM_varimax,T] = rotatefactors(PM_pca_components(:,1:4));

%% Mapping the relaxed solution to MS constrains

M_pca_components = M_pca_components(:,1:4);   %pick the first four components

for i = 1:size(M_pca_components,1)
    vec = M_pca_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_pca_labels(i,:) = vec/max(abs(vec));
    
    vec = M_varmax_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_varimax_labels(i,:) = vec/max(abs(vec));
end

%% Explained variance
M = M';
aprox_pca_norm = M_pca_labels * MS_pca_norm';
aprox_varmax_norm = M_varimax_labels * MS_varmax_norm';

%measured reached residual variance
pca_loadings = diag(M*aprox_pca_norm');
varmax_loadings = diag(M*aprox_varmax_norm');

aprox_pca = aprox_pca_norm.*repmat(pca_loadings,1,19);
aprox_varmax = aprox_varmax_norm.*repmat(varmax_loadings,1,19);
aprox_pca_relaxed = M_pca_components * MS_pca';
aprox_varmax_relaxed = M_varmax_components * MS_varmax';

res_var_pca = trace((M - aprox_pca)*(M - aprox_pca)')/(Nt*(Ns-1));
res_var_varmax = trace((M - aprox_varmax)*(M - aprox_varmax)')/(Nt*(Ns-1));
res_var_pca_relaxed = trace((M - aprox_pca_relaxed)*(M - aprox_pca_relaxed)')/(Nt*(Ns-1));
res_var_varmax_relaxed = trace((M - aprox_varmax_relaxed)*(M - aprox_varmax_relaxed)')/(Nt*(Ns-1));

%theoretical minimum residual variance
min_res_var_varmax = trace(M*M' - (aprox_varmax_norm*M').^2)/(Nt*(Ns-1));
min_res_var_pca = trace(M*M' - (aprox_pca_norm*M').^2)/(Nt*(Ns-1));

%variance of data
var_data = trace(M*M')/(Nt*(Ns-1));

%theoretical explained variance
theo_explained_pca = 1 - min_res_var_pca/var_data;
theo_explained_varimax = 1 - min_res_var_varmax/var_data;

%numerical explained variance
explained_pca = 1 - res_var_pca/var_data;
explained_varmax = 1 - res_var_varmax/var_data;
explained_pca_relaxed = 1 - res_var_pca_relaxed/var_data;
explained_varmax_relaxed = 1 - res_var_varmax_relaxed/var_data;

%print results
disp(['theoretical explained variance by pca: ',num2str(theo_explained_pca)]);
disp(['numerical explained variance by pca: ',num2str(explained_pca)]);
disp(['numerical explained variance by pca relaxed: ',num2str(explained_pca_relaxed)]);
disp(['theoretical explained variance by varimax: ',num2str(theo_explained_varimax)]);
disp(['numerical explained variance by varimax: ',num2str(explained_varmax)]);
disp(['numerical explained variance by varimax relaxed: ',num2str(explained_varmax_relaxed)]);
%% Visualize
figure
suptitle('Varimax Microstates - spatial') %visualize topomaps
for i = 1:size(MS_varmax_norm,2)
    subplot(2,2,i)
    topoplot(MS_varmax_norm(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end

figure
suptitle('PCA Microstates - spatial') %visualize topomaps
for i = 1:size(MS_varmax_norm,2)
    subplot(2,2,i)
    topoplot(MS_pca_norm(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end

figure
suptitle('Varimax Microstates - temporal') %visualize topomaps
for i = 1:size(MS_PM_varimax,2)
    subplot(2,2,i)
    topoplot(MS_PM_varimax(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end

figure
suptitle('PCA Microstates - temporal') %visualize topomaps
for i = 1:size(MS_PM_pca,2)
    subplot(2,2,i)
    topoplot(MS_PM_pca(:,i),'10-20-microst.loc','gridscale',150,'verbose','off');
end
%% Additional visualization

M_pca_sq = M_pca_components.^2;   
M_varimax_sq = M_varmax_components.^2;

