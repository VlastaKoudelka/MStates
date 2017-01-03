%% PCA vs. varimax microstates

close all
clear all
path = 'D:\MicroStates\KubaK\pla008\';
filename = 'PLA_008_EEG_B_r01604-20150223-121157_Edit Channels 2.dat';

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


%% PM approach
PM_pca_components = pca(M_stand');
MS_PM_pca = PM_pca_components(:,1:4);   %here the components are the MS
[MS_PM_varimax,T] = rotatefactors(PM_pca_components(:,1:4));

%% Visualize
figure
suptitle('Varimax Microstates-temporal') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_varmax(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('PCA Microstates-temporal') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_pca(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('Varimax Microstates - spatial') %visualize topomaps
for i = 1:size(MS_PM_varimax,2)
    subplot(2,2,i)
    topoplot(MS_PM_varimax(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('PCA Microstates-spatial') %visualize topomaps
for i = 1:size(MS_PM_pca,2)
    subplot(2,2,i)
    topoplot(MS_PM_pca(:,i),'10-20-microst.loc','gridscale',150);
end

%% Mapping the relaxed solution to MS constrains

M_pca_components = M_pca_components(:,1:4);   %pick the first four components

for i = 1:size(M_pca_components,1)
    vec = M_pca_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_pca_components(i,:) = vec;
    
    vec = M_varimax_components(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_varimax_components(i,:) = vec;
end

%% Explained variance
M_stand = M_stand';
aprox_varimax = M_varimax_components * MS_varmax';
aprox_pca = M_pca_components * MS_pca';

var_varimax = trace((M_stand - aprox_varimax)*(M_stand - aprox_varimax)')/(Nt*(Ns-1))
var_pca = trace((M_stand - aprox_pca)*(M_stand - aprox_pca)')/(Nt*(Ns-1))

var_data = trace(M_stand*M_stand')/(Nt*(Ns-1))

explained_varimax = 1 - var_varimax/var_data
explained_pca = 1 - var_pca/var_data

%% Additional visualization

M_pca_sq = M_pca_components.^2;   
M_varimax_sq = M_varimax_components.^2;




