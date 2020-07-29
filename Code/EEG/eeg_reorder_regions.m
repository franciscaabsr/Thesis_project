brainstorm2free = importdata('/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/destr_brainstorm2freesurfer.txt');
base_dir_eeg = '/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/eeg_connect_destrieux/'; 

for s=1:9
    
    load([base_dir_eeg 'subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_gamma_subj' num2str(s) '-7T_.mat'])
    
    for t=1:length(connEEGgamma)
        
        conn_eeg = connEEGgamma(:,:,t);
        %reorder EEG matrix according to Freesurfer order
        conn_eeg = conn_eeg(brainstorm2free, brainstorm2free);
        connEEGgamma(:,:,t) = conn_eeg;
        
    end
    
    save([base_dir_eeg 'subj0' num2str(s) '-7T/conn_destr_cohi_time_eeg_gamma_correct_subj' num2str(s) '-7T_.mat'], 'connEEGgamma');
      
end