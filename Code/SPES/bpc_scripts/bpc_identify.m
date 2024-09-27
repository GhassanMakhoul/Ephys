function [B_struct,excluded_pairs]=bpc_identify(V,pair_types)
% function [B_struct,excluded_pairs]=bpc_identify(V,pair_types)
% this function identifies basis curves that cluster characteristic shapes
% underlying groups of characteristic shapes (initially implemented with 
% cortico-cortical evoked potentials from ECoG data)
%
% Input variables: 
%
% V: single-trial stimulation-evoked voltage matrix (Dimensions of V are
%    T×K, with K total stimulation events by T total timepoints)
%
% pair_types: structure with subgroups of stimulation pairs, with fields
%    .pair: electrodes associated with the pair, cathode-anode ordering
%    .indices: indices of subgroup in V matrix (out of K total stims)
%
%
% Output variables 
% B_struct: structure with properties of basis curves
%    .curve: shape of basis curve (dimension Tx1), unit normalized
%    .pairs: vector of included stimulation pairs (index into pair_types)
%    .alphas: cell array of alpha values for projections of curve into each stimulation pair subgroup
%    .ep2: cell array of internal projections of residual error for each stimulation trial after application of basis curve
%    .errxproj: cell array of cross projections of residual error for each stimulation trial after application of basis curve
%               (to look for residual structure) 
%
% excluded_pairs: pairs not associated with any basis curve (no residual structure)
% 
% For a full explanation, please see explanatory text & manuscript at: 
% https://purl.stanford.edu/rc201dv0636
% 
% If this code is used in a publication, please cite the manuscript 
% "Basis profile curve identification to understand electrical stimulation 
% effects in human brain networks"
% by KJ Miller, K-R Mueller, and D Hermes
%
% Copyright (C) 2020 Kai J Miller
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% default settings

    % start out number of cluster dimensions at number of pair types (+1)
%         cl_dim=length(pair_types)+1; 
    cl_dim=10;
    
    % number of iterations to re-run NNMF - because can get suboptimal factorization due to non-unique nature of NNMF (starts out with a random matrix)
    %     num_reruns=length(pair_types); 
    num_reruns=20;


%% set of native normalized single stimulation cross-projections
    V0=V./(ones(size(V,1),1)*(sum(V.^2,1).^.5)); % normalize (L2 norm) each trial
    P=V0.'*V; % calculate internal projections

    for k=1:length(pair_types)    
        for l=1:length(pair_types) %2nd index is the normalized one
            if k==l % diagonal 
                % gather all off-diagonal elements from self-submatrix
                a=P(pair_types(k).indices,pair_types(k).indices);
                b=[]; 
                for q=1:(size(a,2)-1), b=[b a(q,(q+1):end)]; end
                for q=2:(size(a,2)), b=[b a(q,1:(q-1))]; end
            else
                b=reshape(P(pair_types(k).indices,pair_types(l).indices),1,[]);
            end
            S{k,l}=b;
            tmat(k,l)=mean(b)/(std(b)/sqrt(length(b))); % calculate t-statistic
        end
    end

    clear k q l a b

%% iteratively decrease inner components of non-negative matrix factorization
    t0=tmat; t0(t0<0)=0; % factorize t-matrix, but first make non-negative
    nnmf_xcorr_score=100; % start off-diagonal penalty score out with high value
    q=0; % index for saving structure
    WH_struct=[]; % saving structure

%     while nnmf_xcorr_score>1 % I chose 1 for sum of off-diagonal elements. Change at your discretion.
    while nnmf_xcorr_score>.5 % Max of off-diagonal elements. Change at your discretion.        
        cl_dim=cl_dim-1; % reduce number of inner dimensions
        q=q+1; % iterate q

        % multiple run-throughs of  nnmf
        tmp_mat=[]; tmp_err=[];
        for k=1:num_reruns
            [tmp_mat(k).W,tmp_mat(k).H,tmp_err(k)] = kjm_nnmf(t0, cl_dim); % without sparseness constraints
        end
        % select factorization with smallest error
        [tmp, k_ind]=min(tmp_err);
        W=tmp_mat(k_ind).W; H=tmp_mat(k_ind).H;

        % score for this decomposition - off diagonal element weights
        % (comment choice in or out)
        
%             % sum of off-diagonals
%             nnmf_xcorr_score=sum(reshape(triu(H*H.',1),1,[])); % sum of nonnegative matrix factorization component interdependencies
%             disp(['Inner dimension ' num2str(cl_dim) ', off diagonal score: ' num2str(nnmf_xcorr_score)])

            % max of off-diagonals
            nnmf_xcorr_score=max(reshape(triu(H*H.',1),1,[])); % sum of nonnegative matrix factorization component interdependencies
            disp(['Inner dimension ' num2str(cl_dim) ', off diagonal maximum element: ' num2str(nnmf_xcorr_score)])
        
        % save matrices and scores for plotting later
        WH_struct(q).W=W;
        WH_struct(q).H=H;
        WH_struct(q).nnmf_xcorr_score=nnmf_xcorr_score;
    end

%% winner take all on columns of H and then threshold by 1/(2*sqrt(length(pair_types))) -- since all equal would be 1/sqrt(length(pair_types))
    H0=0*H;
    for k=1:length(pair_types)
        [tmp,k_ind]=max(H(:,k));
        H0(k_ind,k)=tmp;
    end

    % Here is where we can change threshold
    H0=H0>1/(2*sqrt(length(pair_types)));
%     H0=H0>1/(sqrt(length(pair_types)));
%       H0=H0>1/(3*sqrt(length(pair_types)));
  
%% now can find clusters of stimulation pair groups from those that aren't dramatically distinguishable from one another (after BF correction)
    for q=1:size(H0,1)    
       
        cl_pg=find(H0(q,:)); % cluster pair groups
        cl_inds=[]; % cluster indices (e.g. concatenated trials from cluster pair groups)
        for k=1:length(cl_pg)
            cl_inds=[cl_inds; pair_types(cl_pg(k)).indices];
        end
        V_tmp=V(:,cl_inds);

        % use the "kernel trick" to estimate eigenvectors of this cluster of pair groups
        %  keep 1st eigenvector for basis
        [E,tmp]=kt_pca(V_tmp);
        B_tmp=E(:,1); % basis vector is 1st PC
        if mean(B_tmp.'*V_tmp)<0 % invert if it projects into its own group negatively
            B_tmp=-B_tmp;
        end

        B(q,:)=B_tmp; % matrix of basis vectors
        B_struct(q).curve=B_tmp;
        B_struct(q).pairs=cl_pg;
    end
    
    % pairs not represented by any basis
    excluded_pairs=find(1-sum(H0,1));
    
%% calculate statistics for each basis curve, eliminating those where there is no significant representation
    for bb=1:size(B,1) % cycle through basis curves

        al=B(bb,:)*V; % alpha coefficient weights for basis curve bb into V
        ep=V-B(bb,:).'*al; % residual epsilon (error timeseries) for basis bb after alpha*B coefficient fit
        errxproj=ep.'*ep; % calculate all projections of error
        V_selfproj=V.'*V; % power in each trial

    %    
        for n=1:length(B_struct(bb).pairs)   % cycle through pair types represented by this basis curve
            tmp_inds=pair_types(B_struct(bb).pairs(n)).indices; %indices for this pair type

            B_struct(bb).alphas{n}=al(tmp_inds); % alpha coefficient weights for basis curve bb into V

            a=errxproj(tmp_inds,tmp_inds); % self-submatrix of error projections   

            B_struct(bb).ep2{n}=diag(a).'; % sum-squared error
            B_struct(bb).V2{n}=diag(V_selfproj(tmp_inds,tmp_inds)).'; % sum-squared individual trials
            

            % gather all off-diagonal elements from self-submatrix
            b=[]; 
            for q=1:(size(a,2)-1), b=[b a(q,(q+1):end)]; end
            for q=2:(size(a,2)), b=[b a(q,1:(q-1))]; end        
            B_struct(bb).errxproj{n}=b; % systematic residual structure within a stim pair group for a given basis will be given by set of native normalized internal cross-projections
            
            
            clear q a b tmp_inds
        end
        clear n al ep errxproj
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [W,H, newerr] = kjm_nnmf(V,rdim)
% kjm_nnmf - non-negative matrix factorization 
% 
% SYNTAX:
% [W,H] = kjm_nmfsc(V, rdim);
%
% INPUTS:
% V          - data matrix (N x T)
% rdim       - (M) number of 'reduced dimension' components (inner dimension of factorization)
% 
% OUTPUTS:
% W           - first non-negative factor in V~WH decomposition (N x rdim)
% H           - second non-negative factor in V~WH decomposition (rdim x T)
% newerr      - latest estimate of error of in estimation
%
% file is based upon elements of file nmfsc by Patrik Hoyer at https://github.com/aludnam/MATLAB/tree/master/nmfpack
% ref is Hoyer, 2004, 'Non-negative matrix factorization with sparseness constraints' Journal of Machine Learning Research  5:1457-1469, 2004.
% kjm 5/2020


%% defaults and data check

    conv_thresh=1e-5; %convergence threshold

    % Check that we have non-negative data
    if min(V(:))<0, error('Negative values in data!'); end

    % Globally rescale data to avoid potential overflow/underflow
    V = V/max(V(:));
    
%% Create initial random positive matrices
    W = abs(randn(size(V,1),rdim)); 
    H = abs(randn(rdim,size(V,2)));
    
    % normalizes so each row of H has unit L2 norm
    H = H./(sqrt(sum(H.^2,2))*ones(1,size(V,2)));

%% initialize iterative process
    % Calculate initial error (E)
    errhistory = 0.5*sum(sum((V-W*H).^2));
    loop_continue=1;
    loop_num=0;
    
%% iterative convergence
while loop_continue==1

    loop_num=loop_num+1;
        
    %% ----- Update H ---------------------------------------
        % Update using standard NMF multiplicative update rule, (element-wise multiplication of ratios of positive numbers ensures positivity)
        H = H.*(W.'*V)./(W.'*W*H + 1e-9); %adds 1e-9 so doesn't explode if divide by zero

        % Renormalize so rows of H have unit norm
        Hnorms = sqrt(sum(H'.^2));
        H = H./(Hnorms'*ones(1,size(V,2)));
        W = W.*(ones(size(V,1),1)*Hnorms); % pushes scaling from rows of H onto columns of W
	    
    %% ----- Update W ---------------------------------------
        % Update using standard NNMF multiplicative update rule 
        % (element-wise multiplication of ratios of positive numbers ensures positivity)
        W = W.*(V*H.')./(W*H*H.' + 1e-9);	
	        
    %% ----- Calculate error (E) ---------------------------------------
        newerr = norm(V-W*H,'fro')^2;
        errhistory = [errhistory newerr];

        steperrratio=abs(diff(errhistory([-1 0]+end)))/newerr;
        if steperrratio<conv_thresh, 
%             disp(['ratio of change in error to error below threshold of ' num2str(conv_thresh)]);
            loop_continue=0;
        end
%     
%         if or(mod(loop_num,1000)==1,loop_continue==0)
%             disp(['Iteration: ' num2str(loop_num)...
%                 ', error term = ' num2str(newerr)...
%                 ', error step ratio = ' num2str(steperrratio)])
%         end
end

[tmp,pp]=sort(sum(W.^2,1),'descend');
W=W(:,pp);H=H(pp,:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [E,S]=kt_pca(X)
% function [E,S]=kt_pca(X)
% This is an implementation of the linear kernel PCA method ("kernel trick")
% described in "Kernel PCA Pattern Reconstruction via Approximate Pre-Images"
% by Schölkopf et al, ICANN, 1998, pp 147-15
% See also course lectures by K.R. Muller (TU Berlin)
% 
% Inputs:
% X(T,N) - Matrix of data in. Only need this trick if T>>N
% 
% Outputs:
% E(T,N) - Columns of E are estimated Eigenvectors of X, in descending order
% S(1,N) - Corresponding estimated Eigenvalues of X, in descending order
% 
% kjm, june 2020

%% use the "kernel trick" to estimate eigenvectors of this cluster of pair groups
    [F,S2]=eig(X.'*X); % eigenvector decomposition of (covariance of transpose) - may want to replace this with cov function so mean is subtracted off
    [S2,v_inds]=sort(sort(sum(S2)),'descend'); F=F(:,v_inds); %reshape properly    

    S=S2.^.5; % estimated eigenvalues of both X.'*X and X*X.'
    
    %     ES=X*F.';
    ES=X*F; % kernel trick
    E=ES./(ones(size(X,1),1)*S); % divide through to obtain unit-normalized eigenvectors



