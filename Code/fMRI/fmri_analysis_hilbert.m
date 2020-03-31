function fmri_analysis_hilbert(mode,state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analysis of fMRI data using Leading Eigenvector Dynamics Analysis (LEiDA)
% Code adapted from Joana Cabral
%
% This function processes, clusters and analyses BOLD data using LEiDA
%
% Note: Step 4 (display of all results) can be run independently once
%  fmri_results_Hilbert.mat has been saved:
%  Choose mode = 'run' or 'read' to run LEiDA, or read after results are saved
%  In the 'run' mode, the clusters and the probabilites, lifetimes and
%  transition matrices associated to each dFC state will be computed for
%  all states in a range from 3 to 15 states (can be changed)
%  In the 'read' mode, specify the number of states as a second argument, so that
%  only the results from that number of states are displayed
%
% 1 - Read the BOLD data from the folders and computes the BOLD phases
%   - Calculate the instantaneous BOLD synchronization matrix
%   - Compute the Leading Eigenvector at each frame from all fMRI scans
% 2 - Cluster the Leading Eigenvectors into recurrent Functional Networks
% 3 - Compute the probability, lifetimes and transition matrices of each FN
%   - Saves the Eigenvectors, Clusters and statistics into fmri_results_Hilbert.mat
%
% 4 - Plots FC states for each clustering solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(mode,'run')
    %% Initialization
    
    % Data parameters:
    N_subjects = 9;
    N_areas    = 90;
    fmri_TR    = 1;
    Tmax       = 470;
    
    % Folder where data is saved.
    Directory  = '/home/francisca/Documents/Tese/Implementation/Datasets/BOLDmatrices/';
    
    %% Obtain phase coherence matrices using Hilbert transform
    
    % Bandpass filter settings
    fnq = 1/(2*fmri_TR);          % Nyquist frequency
    flp = .01;                    % lowpass frequency of filter (Hz)
    fhi = .1;                     % highest BOLD frequency considered (Hz)
    
    Wn = [flp/fnq fhi/fnq];       % butterworth bandpass non-dimensional frequency
    k = 2;                        % 2nd order butterworth filter
    [bfilt,afilt] = butter(k,Wn); % construct the filter
    clear fnq flp fhi Wn k
    
    % Preallocate variables to save FC patterns and associated information
    
    V1_all   = zeros(Tmax*N_subjects,N_areas); % All leading eigenvectors
    t_all    = 0; % Index of time (starts at 0 and will be updated until N_subj)
    Time_all = zeros(Tmax*N_subjects,1); % Vector to indicate the subject number at each TR
    Var_Eig  = zeros(N_subjects,Tmax); % Saves the variance explained by the leading eigenvectors
    iFC_all  = zeros(N_areas,N_areas,Tmax,N_subjects); % 4D matrix to store PC for each TR and subject
    
    for s = 1:N_subjects
        
        % USER: Adapt to load here the BOLD matrix (NxT) from each subject
        load([Directory 'BOLD' num2str(s) '.mat'])
        BOLD = tc_aal; % N_areas x Tmax
        
        % Filter the BOLD & get the BOLD phase using the Hilbert transform
        
        Phase_BOLD = zeros(N_areas,Tmax);
        
        for seed = 1:N_areas
            BOLD(seed,:) = BOLD(seed,:)-mean(BOLD(seed,:));
            signal_filt = filtfilt(bfilt,afilt,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        for t=1:Tmax
            
            %Calculate the Instantaneous FC (BOLD Phase Synchrony)
            iFC = zeros(N_areas);
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            iFC_all(:,:,t,s) = iFC;
            
            % Get the leading eigenvector
            
            [eVec,eigVal]=eigs(iFC);
            eVal=diag(eigVal);
            [val1, i_vec_1] = max(eVal);
            V1=eVec(:,i_vec_1);
            if sum(V1)<0
                V1=-V1;
            end
            
            % Compute the variance explained by the leading eigenvector
            Var_Eig(s,t)=val1/sum(eVal);
            
            % Save V1 from all frames in all fMRI sessions
            t_all=t_all+1; % Update time
            V1_all(t_all,:)=V1;
            Time_all(t_all)=s;
        end
        
    end
    
    %% Cluster the Leading Eigenvectors
    
    disp('Starting clustering')
    
    maxk=15;
    mink=3;
    rangeK=mink:maxk;
    
    % Set the parameters for Kmeans clustering
    Kmeans_results=cell(size(rangeK));
    
    for k=1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' FC states'])
        [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance','sqeuclidean','Replicates',1000,'MaxIter',300,'Display','final','Options',statset('UseParallel',1));
        [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
        [~,idx_sort]=sort(ind_sort,'ascend');
        Kmeans_results{k}.IDX=idx_sort(IDX);    % Cluster time course - numeric collumn vectors
        Kmeans_results{k}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
        Kmeans_results{k}.SUMD=SUMD(ind_sort);  % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D(:,ind_sort);       % Distance from each point to every centroid
    end
    
    %% Smoothing dFC labels
    
    % Smoothing is done in order to impose a temporal restriction (as k-means
    % does not do it)
    
    clear par
    par.b = 10; % window size chose based on the higher limit of band pass filter (0.1 Hz -> 10 s)
    par.lambda = 0.5;
    dFC_labels = zeros(length(rangeK),Tmax*N_subjects);
    for i=1:length(rangeK)
        for s = 1:N_subjects
            dFC_labels(i,(s-1)*470+1:s*470)= SmoothLabelsKmeans(V1_all((s-1)*470+1:s*470,:)',Kmeans_results{i}.C,par,1/fmri_TR);
        end
    end
    
    
    %% mean lifetime for all rangeK
    
    lifetime_allk = cell(1,length(rangeK));
    for i= 1:length(rangeK)
        lifetime = zeros(1,rangeK(i));
        stderr_lifetime = zeros(1,rangeK(i));
        for k=1:rangeK(i)
            lifetime_subj=[];
            for s=1:N_subjects
                Ctime_bin=dFC_labels(i,(s-1)*Tmax+1:s*Tmax)==k;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % Discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                    lifetime_subj=[lifetime_subj C_Durations];
                end
            end
            lifetime(k) = mean(lifetime_subj);
            stderr_lifetime(k) = std(lifetime_subj)/sqrt(length(lifetime_subj));
        end
        lifetime_allk{i}.mean = lifetime;
        lifetime_allk{i}.err = stderr_lifetime;
    end
    
    %% Fractional Occupancy for all rangeK
    
    FO_all = cell(1,length(rangeK));
    for i=1:length(rangeK)
        FO=zeros(rangeK(i),1);
        for k=1:rangeK(i)
            FO(k)=mean(dFC_labels(i,:)==k);
        end
        FO_all{i}=FO;
    end
    
    %% Transition matrix for all rangeK
    
    Trans_mat_all = cell(1,length(rangeK));
    for i = 1:length(rangeK)
        trans_mat=zeros(rangeK(i));
        idx = dFC_labels(i,:);
        switch_state=diff(idx);
        for j=1:length(idx)-1
            if switch_state(j) ~= 0
                trans_mat(idx(j),idx(j+1)) = trans_mat(idx(j),idx(j+1)) + 1;
            end
        end
        soma=sum(trans_mat,2);
        for k=1:rangeK(i)
            if soma(k) ~= 0
                trans_mat(k,:)=trans_mat(k,:)/soma(k);
            end
        end
        Trans_mat_all{i} = trans_mat;
    end
    
    save('fmri_results_Hilbert.mat','V1_all','Var_Eig','Kmeans_results','rangeK','N_subjects','dFC_labels','lifetime_allk','FO_all','Trans_mat_all','N_areas','iFC_all')
    
elseif strcmp(mode,'read')
    
    load('fmri_results_Hilbert.mat')
    
    state_idx = state-min(rangeK)+1;
    
    %% Display mean lifetime result
    
    figure
    bar(lifetime_allk{state_idx}.mean)
    ylabel('Mean lifetime (TR)','fontsize',18)
    xlabel('State','fontsize',18)
    title('Mean lifetime','fontsize',20)
    set(gca,'fontsize',16)
    
    hold on
    
    errorbar(1:state, lifetime_allk{state_idx}.mean, lifetime_allk{state_idx}.err,'LineStyle','none','Color','k');
    
    hold off
    %% Display Fractional Occupancy result
    
    figure
    bar(FO_all{state_idx})
    xlabel('State','fontsize',18)
    ylabel('Fractional Occupancy','fontsize',18)
    title('Fractional Occupancy','fontsize',20)
    set(gca,'fontsize',16)
    %% Display transition matrix result
    
    figure
    colormap;
    trans_mat = Trans_mat_all{state_idx};
    imagesc(trans_mat)
    colorbar
    title('Transition Matrix','fontsize',20)
    ylabel('From State','fontsize',18)
    xlabel('To State','fontsize',18)
    set(gca,'xtick',1:14,'ytick',1:14,'fontsize',16)
    thresh = find(trans_mat > 0.15);
    textStrings = num2str(trans_mat(thresh), '%0.2f');       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:state);  % Create x and y coordinates for the strings
    hStrings = text(x(thresh), y(thresh), textStrings(:),'HorizontalAlignment', 'center');
    
    %% Display nodes in cortex
    
    V = Kmeans_results{state_idx}.C;
    %Order=[1:2:N_areas N_areas:-2:2];
    
    figure('defaultaxesfontsize',12)
    colormap(jet)
    for c=1:state
    
        subplot(2,state,c)
        plot_nodes_in_cortex(V(c,:))
        title({['State #' num2str(c)]},'fontsize',18)
    
        subplot(2,state,c+state)
        Vc=V(c,:);
        Vc=Vc/max(abs(Vc));
        FC_V=Vc'*Vc;
    
        li=max(abs(FC_V(:)));
        imagesc(FC_V,[-li li])
        axis square
        if c ==1
            ylabel('Brain area #','fontsize',12)
            xlabel('Brain area #','fontsize',12)
        end
        set(gca,'xtick',0:45:90,'ytick',0:45:90)
    
    end
    
    for c=1:state
        
        figure
        colormap(jet)
        subplot(2,1,1)
        plot_nodes_in_cortex(V(c,:))
        title({['State #' num2str(c)]},'fontsize',18)
        
        subplot(2,1,2)
        Vc=V(c,:);
        Vc=Vc/max(abs(Vc));
        FC_V=Vc'*Vc;
        
        li=max(abs(FC_V(:)));
        imagesc(FC_V,[-li li])
        axis square
        colorbar
        ylabel('Brain area #','fontsize',10)
        xlabel('Brain area #','fontsize',10)
        set(gca,'xtick',0:45:90,'ytick',0:45:90,'fontsize',14)
        
    end
    
    
    load AAL_labels.mat
    
    figure
    for c = 1:state
        Vc=V(c,:);
        Vc=Vc/max(abs(Vc));
        subplot(1,state,c)
        hold on
        barh(Vc.*(Vc<=0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
        barh(Vc.*(Vc>0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
        ylim([0 91])
        xlim([-1 1])
        grid on
        set(gca,'YTick',1:N_areas,'Fontsize',6,'YDir','reverse')
        if c == 1
            set(gca,'YTickLabel',[num2str([1:90]') blanks(90)' label90])
        else
            set(gca,'YTickLabel',1:90)
        end
    end
end
end