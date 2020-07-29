for s=1:9
    
    load(['eeg_connect_destrieux/subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_gamma_correct_subj' num2str(s) '-7T_.mat'])
    
    conn_av = zeros(148,148);
    %aux_fmri = zeros(68,68);

    conn_av = mean(connEEGgamma,3);
    %for t = 1:length(connEEGgamma)
    %    conn_av = conn_av + connEEGgamma(:,:,t);

        %[ii,jj] = find(connEEGgamma(:,:,t)); %get index of nonzero values on connectivity matrix
        %aux = zeros(68,68);
        %idx = sub2ind(size(aux), ii, jj);
        %if ~isempty(idx)
        %    aux(idx) = 1;
        %end
        %aux_fmri = aux_fmri + aux; %add count for the nonzero values
    %end

    %[ii,jj] = find(~aux_fmri);
    %idx = sub2ind(size(aux_fmri), ii, jj);
    %aux_fmri(idx) = 1; % so as not to divide by zero

    %conn_av = rdivide(conn_av,aux_fmri);
    %conn_av = conn_av/length(connEEGgamma);
    %conn_av = conn_av - diag(diag(conn_av));

    fig = figure;
    set(gcf,'Units','inches', 'Position',[0 0 6 4])

    %clim = [ -1 1 ];
    %im = imagesc(conn_av, clim);
    im = imagesc(conn_av);
    colormap(jet);
    %colormap(flipud(hot(10)));

    h = colorbar('eastoutside');
    xlabel(h, 'h', 'FontSize', 14);

    % Title, axis
    title(['Average connectivity matrix' num2str(s) 'EEG gamma'], 'FontSize', 14);
    set(gca, 'Ticklength', [0 0])
    grid off
    box off
    
    name = sprintf('/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/static_connect_matrices_destrieux/subj0%s/average_eeg_gamma_correct.png', num2str(s));
    saveas(fig,name);

end