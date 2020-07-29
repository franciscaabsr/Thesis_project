%regions = 68;
regions = 148;
%% STATIC CONNECTOME FMRI

mrtx = zeros(regions);
count = 1;
for r1 = 1:regions-1
    for r2 = r1+1:regions
        mrtx(r1,r2) = subj(9).sess.fMRI(count);
        mrtx(r2,r1) = mrtx(r1,r2);
        count = count + 1;
    end
end

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])

clim = [ -1 1 ];
%im = imagesc(conn(19:86,19:86), clim);
im = imagesc(mrtx,clim);

colormap(jet);

h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);

% Title, axis
title('Static connectome fMRI subj09 - Jonathan ', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

%% STATIC CONNECTOME EEG 

for s=1:9
    mrtx = zeros(regions);
    count = 1;
    for r1 = 1:regions-1
        for r2 = r1+1:regions
            mrtx(r1,r2) = subj(s).sess.fMRI(count);
            mrtx(r2,r1) = mrtx(r1,r2);
            count = count + 1;
        end
    end

    fig = figure;
    set(gcf,'Units','inches', 'Position',[0 0 6 4])

    %clim = [ 0 1 ];
    %im = imagesc(mrtx, clim);
    im = imagesc(mrtx);
    colormap(jet);

    h = colorbar('eastoutside');
    xlabel(h, 'h', 'FontSize', 14);

    % Title, axis
    title('Static connectome fMRI - Jonathan ', 'FontSize', 14);
    set(gca, 'Ticklength', [0 0])
    grid off
    box off
    
    name = sprintf('/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/subj0%s/static_fmri.png', num2str(s));
    saveas(fig,name);
    
end


%% STATIC COMPARISON - concatenating all subjects

conn_av_fmri = zeros(regions,regions);
conn_av_broad = zeros(regions,regions);
conn_av_alpha = zeros(regions,regions);
conn_av_beta = zeros(regions,regions);
conn_av_delta = zeros(regions,regions);
conn_av_gamma = zeros(regions,regions);
conn_av_theta = zeros(regions,regions);

num_time_points = 0;

% to store number of times non zero
aux_fmri = zeros(regions,regions);
aux_broad = zeros(regions,regions);
aux_alpha = zeros(regions,regions);
aux_beta = zeros(regions,regions);
aux_delta = zeros(regions,regions);
aux_gamma = zeros(regions,regions);
aux_theta = zeros(regions,regions);

%get average connectivity matrices, corresponding to the static connectome
for s = 1:9
    
    load(['fmri_connect_destrieux/subj0' num2str(s) '-7T/conn_destrx_phase_coh_time_fmri.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_broad_correct_subj' num2str(s) '-7T_.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_alpha_correct_subj' num2str(s) '-7T_.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_beta_correct_subj' num2str(s) '-7T_.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_delta_correct_subj' num2str(s) '-7T_.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_gamma_correct_subj' num2str(s) '-7T_.mat'])
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_theta_correct_subj' num2str(s) '-7T_.mat'])
    
    num_time_points = num_time_points + length(connFMRI);
    
    conn_av_fmri = conn_av_fmri + mean(connFMRI,3);
    conn_av_broad = conn_av_broad + mean(connEEGbroad,3);
    conn_av_alpha = conn_av_alpha + mean(connEEGalpha,3);
    conn_av_beta = conn_av_beta + mean(connEEGbeta,3);
    conn_av_delta = conn_av_delta + mean(connEEGdelta,3);
    conn_av_gamma = conn_av_gamma + mean(connEEGgamma,3);
    conn_av_theta = conn_av_theta + mean(connEEGtheta,3);

%     for t = 1:length(connFMRI)
%         
%         conn_av_fmri = conn_av_fmri + connFMRI(:,:,t);
%         [ii,jj] = find(connFMRI(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_fmri = aux_fmri + aux; %add count for the nonzero values
%         
%         conn_av_broad = conn_av_broad + connEEGbroad(:,:,t);
%         [ii,jj] = find(connEEGbroad(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_broad = aux_broad + aux; %add count for the nonzero values
%         
%         conn_av_alpha = conn_av_alpha + connEEGalpha(:,:,t);
%         [ii,jj] = find(connEEGalpha(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_alpha = aux_alpha + aux; %add count for the nonzero values
%         
%         conn_av_beta = conn_av_beta + connEEGbeta(:,:,t);
%         [ii,jj] = find(connEEGbeta(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_beta = aux_beta + aux; %add count for the nonzero values
%         
%         conn_av_delta = conn_av_delta + connEEGdelta(:,:,t);
%         [ii,jj] = find(connEEGdelta(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_delta = aux_delta + aux; %add count for the nonzero values
%         
%         conn_av_gamma = conn_av_gamma + connEEGgamma(:,:,t);
%         [ii,jj] = find(connEEGgamma(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_gamma = aux_gamma + aux; %add count for the nonzero values
%         
%         conn_av_theta = conn_av_theta + connEEGtheta(:,:,t);
%         [ii,jj] = find(connEEGtheta(:,:,t)); %get index of nonzero values on connectivity matrix
%         aux = zeros(regions,regions);
%         idx = sub2ind(size(aux), ii, jj);
%         aux(idx) = 1;
%         aux_theta = aux_theta + aux; %add count for the nonzero values
%         
%     end
end

% [ii,jj] = find(~aux_fmri);
% idx = sub2ind(size(aux_fmri), ii, jj);
% aux_fmri(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_broad);
% idx = sub2ind(size(aux_broad), ii, jj);
% aux_broad(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_alpha);
% idx = sub2ind(size(aux_alpha), ii, jj);
% aux_alpha(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_beta);
% idx = sub2ind(size(aux_beta), ii, jj);
% aux_beta(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_delta);
% idx = sub2ind(size(aux_delta), ii, jj);
% aux_delta(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_gamma);
% idx = sub2ind(size(aux_gamma), ii, jj);
% aux_gamma(idx) = 1; % so as not to divide by zero
% 
% [ii,jj] = find(~aux_theta);
% idx = sub2ind(size(aux_theta), ii, jj);
% aux_theta(idx) = 1; % so as not to divide by zero
% 
% conn_av_fmri = rdivide(conn_av_fmri,aux_fmri);
% conn_av_broad = rdivide(conn_av_broad,aux_broad);
% conn_av_theta = rdivide(conn_av_theta,aux_theta);
% conn_av_theta = rdivide(conn_av_theta,aux_theta);
% conn_av_theta = rdivide(conn_av_theta,aux_theta);
% conn_av_theta = rdivide(conn_av_theta,aux_theta);
% conn_av_theta = rdivide(conn_av_theta,aux_theta);

conn_av_fmri = conn_av_fmri/9;
conn_av_broad = conn_av_broad/9;
conn_av_alpha = conn_av_alpha/9;
conn_av_beta = conn_av_beta/9;
conn_av_delta = conn_av_delta/9;
conn_av_gamma = conn_av_gamma/9;
conn_av_theta = conn_av_theta/9;

%set diagonal to zeros as the EEG connectivity matrices
conn_av_fmri = conn_av_fmri - diag(diag(conn_av_fmri));

%Plot colourful matrices 
fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_fmri);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix fMRI', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_broad);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG broad', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_alpha);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG alpha', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_beta);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG beta', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_delta);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG delta', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_gamma);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG gamma', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

fig = figure;
set(gcf,'Units','inches', 'Position',[0 0 6 4])
im = imagesc(conn_av_theta);
colormap(jet);
h = colorbar('eastoutside');
xlabel(h, 'h', 'FontSize', 14);
title('Average connectivity matrix EEG theta', 'FontSize', 14);
set(gca, 'Ticklength', [0 0])
grid off
box off

% connectivity between fMRI and EEG

A = triu(conn_av_fmri);
vec_fmri = nonzeros(A.');

%comparison static connectomes - fMRI vs EEG broad
A = triu(conn_av_broad);
vec_broad = nonzeros(A.');

r = corrcoef(vec_fmri,vec_broad);
fig = figure;
im = scatter(vec_fmri,vec_broad);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG broad connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_broad_correct_dstrx.png')

%comparison static connectomes - fMRI vs EEG alpha
A = triu(conn_av_alpha);
vec_alpha = nonzeros(A.');

r = corrcoef(vec_fmri,vec_alpha);
fig = figure;
im = scatter(vec_fmri,vec_alpha);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG alpha connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_alpha_correct_dstrx.png')

%comparison static connectomes - fMRI vs EEG beta
A = triu(conn_av_beta);
vec_beta = nonzeros(A.');

r = corrcoef(vec_fmri,vec_beta);
fig = figure;
im = scatter(vec_fmri,vec_beta);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG beta connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_beta_correct_dstrx.png')

%comparison static connectomes - fMRI vs EEG delta
A = triu(conn_av_delta);
vec_delta = nonzeros(A.');

r = corrcoef(vec_fmri,vec_delta);
fig = figure;
im = scatter(vec_fmri,vec_delta);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG delta connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_delta_correct_dstrx.png')

%comparison static connectomes - fMRI vs EEG gamma 
A = triu(conn_av_gamma);
vec_gamma = nonzeros(A.');

r = corrcoef(vec_fmri,vec_gamma);
fig = figure;
im = scatter(vec_fmri,vec_gamma);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG gamma connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_gamma_correct_dstrx.png')

%comparison static connectomes - fMRI vs EEG theta
A = triu(conn_av_theta);
vec_theta = nonzeros(A.');

r = corrcoef(vec_fmri,vec_theta);
fig = figure;r
im = scatter(vec_fmri,vec_theta);
xlabel('fMRI connectivity', 'FontSize', 10);
ylabel('EEG theta connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.3f',r(2)));

saveas(fig,'/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/comparison_fMRI_EEG_frequency_bands/fmri_vs_theta_correct_dstrx.png')


%% test

conn_test = zeros(10878,1);

for s=1:9
    
    conn_test = conn_test + subj(s).sess.theta;
    
end

conn_test = conn_test/9;

figure;
r = corrcoef(vec_theta,conn_test);
im = scatter(vec_theta,conn_test);
xlabel('EEG theta connectivity', 'FontSize', 10);
ylabel('EEG theta Jonathan connectivity', 'FontSize', 10);
text(0.6,0.1,sprintf('r = %0.2f',r(2)));

