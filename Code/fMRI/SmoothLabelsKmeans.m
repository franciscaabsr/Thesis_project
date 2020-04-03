%SmoothLabelsKmeans() Smoothes k-means labels
%
% Usage:
%   >> smoothedLabels = SmoothLabelsKmeans(V,M,par, srate)
%
% Where: - V is the input vector on k-means algorithm (Ne x Nt)
%        - M is the centroid vector of each cluster(Nu x Ne)
%        - par.b is the window size
%        - par.lambda is the non-smoothness penalty factor
%          (should be between 0 and 1)
%        - srate is the sampling rate
%        
% Adapted from SmoothLabels from Microstate EEGlab toolbox

function smoothedLabels = SmoothLabelsKmeans(V,M,par, srate)
    
    if ~isfield(par,'b');   par.b = 0;  end
    [Ns,Nt] = size(V);
    Nu = size(M,1);
    deltaT = 1 / srate;
    WinPts = round(par.b / deltaT);
    Vvar = sum(V.*V,1);
    rmat = repmat((1:Nu)',1,Nt);
        
    %Trivial case with one map 
    if(Nu == 1);    smoothedLabels = ones(1,Nt);    return; end
 
    % Computation of squared euclidean distance between each eigenvector and each
    % cluster centroid
    
    fit = zeros(Nu,Nt);
    for i=1:Nu
        for j=1:Nt
            fit(i,j) = norm(M(i,:)-V(:,j)');
        end
    end
       
    [~,smoothedLabels] = min(fit);
    
    % No smoothing
    if WinPts == 0
        w = zeros(Nu,Nt);
        w(rmat == repmat(smoothedLabels,Nu,1)) = 1;  
        return; 
    end
    
    % Some helpful stuff

    crit = 10e-6;
    S0 = 0;
    
    
    % Line 3 and 4
    w = zeros(Nu,Nt);
    w(rmat == repmat(smoothedLabels,Nu,1)) = 1;
    
    e = sum(sum(fit.*w))/(Nt * (Ns - 1));

    DoLoop = true;

    while DoLoop == true
        % Line 5a
        Nb = conv2(w,ones(1,WinPts+1),'same');
        
        % Line 5b
         x = (fit) / (2* e * (Ns-1)) - par.lambda * Nb;
        [~,dlt] = min(x,[],1);
        
        % Line 6
        smoothedLabels = dlt;

        % Line 7
        w = zeros(Nu,Nt);
        w(rmat == repmat(smoothedLabels,Nu,1)) = 1;
        Su = sum(sum(fit.*w)) / (Nt * (Ns - 1));
        
        % Line 8
        if abs(Su - S0) <= crit * Su
            DoLoop = false;
        end
    
        S0 = Su;
    end
    smoothedLabels(Vvar == 0) = 0;
end
