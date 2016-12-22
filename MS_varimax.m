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
M_pca = pca(M_stand);       %initial factors from PCA
MS_pca = M_stand * M_pca(:,1:4); %project original data into PCA space
[M_varimax,T] = rotatefactors(M_pca(:,1:4)); %rotate the first four (varimax is default)
MS_varmax = M_stand * M_varimax;                  %project original data into rotated space

figure
suptitle('Varimax Microstates') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_varmax(:,i),'10-20-microst.loc','gridscale',150);
end

figure
suptitle('PCA Microstates') %visualize topomaps
for i = 1:size(MS_varmax,2)
    subplot(2,2,i)
    topoplot(MS_pca(:,i),'10-20-microst.loc','gridscale',150);
end

%% Mapping the relaxed solution to MS constrains

M_pca = M_pca(:,1:4);   %pick the first four components

for i = 1:size(M_pca,1)
    vec = M_pca(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_pca(i,:) = vec;
    
    vec = M_varimax(i,:);
    vec(abs(vec)<max(abs(vec))) = 0;
    M_varimax(i,:) = vec;
end

%% Explained variance
M_stand = M_stand';
aprox_varimax = M_varimax * MS_varmax';
aprox_pca = M_pca * MS_pca';

var_varimax = trace((M_stand - aprox_varimax)*(M_stand - aprox_varimax)')/(Nt*(Ns-1))
var_pca = trace((M_stand - aprox_pca)*(M_stand - aprox_pca)')/(Nt*(Ns-1))

var_data = trace(M_stand*M_stand')/(Nt*(Ns-1))

explained_varimax = 1 - var_varimax/var_data
explained_pca = 1 - var_pca/var_data

%% Additional visualization

M_pca_sq = M_pca.^2;   
M_varimax_sq = M_varimax.^2;




