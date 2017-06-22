function outp = CDAA(tfp, tffp, brdr, n_clst, cl_cntrs, goi, rt, thrs)
%CDAA predicts regulatory connections between genes based on time course
%transcriptomic data.
%
%   The function can be run with no parameter at the beginning 
%       (Ex.:>>CDAA). Messages and outputs produced by the function will
%       lead to the selection of input values to provide.
%
%   outp - variable that contains either the next input value to the
%       function during the variable selection process or the predicted 
%       collection of regulators/targets when all variables are set.
%
%   tfp - path to the file with transcriptome time course data. The file
%       should be in .csv format where rows correspond to differentially
%       expressed genes and columns to sampling time points. The first row
%       contains values of sampling time points, first column contains gene
%       identifiers (agi's). 
%       Example of .csv file opened in MS Excell might look as follows:
%       |         |   1   |   3   | ... |  24  |
%       |At1g23456| 1.223 | 2.334 | ... | 0.44 |
%       |   ...   |  ...  |  ...  | ... | .... |
%       |At3g45678| 0.523 | 3.224 | ... | 2.14 | 
%
%   tffp - path to the file with agi's of genes that are known to be
%       transcription factors. The file should be in .csv and contain
%       genes' agi's in the first column with no empty lines in between. 
%       This variable can be set as [] if no such file exist. In that case
%       all genes would be assumed to be TFs.
%   
%   brdr - an array containing the border between Initiation and responce
%       stages value 'b' and the number of time points in the Primary 
%       response stage 'c' in the formst [b,c].
%
%   n_clst - an array containing the selected number of clusters at each
%       except for the last one stages. this variable is only relevant
%       during the cluster number selection process. After that numbers of
%       clusters at each state are defined by the number of cluster
%       centroids in 'cl_cntrs'.
%
%   cl_cntrs - a cell array with cluster centroids obtained during cluster
%       number selection process. Example for 3 stages with 2 clusters: 
%       cl_cntrs = {[1.2, 0.5; 0.2, -0.3], [0, 0, 0.7; 0, 0, -0.2]}
%
%   goi - agi of the gene of interest.
%
%   rt - the type of regulatory connections to discover for the gene of
%       interest: set 'regulators' to obtain possible regulators 
%       for the gene of interest or 'targets' to obtain possible targets
%
%   thrs - an array of threshold values. First value corresponds to a
%       threshold for dissimilarity values, the rest - thresholds for
%       denoising process. If not set, the default value id [0.4 0.2 0.2].

    % Initializing output
    outp = 0;
    
    % First start with no parameters provided
    if nargin < 1
        fprintf('\nSet the path to the transcriptomic data file\n');
        fprintf(' as the first parameter\n');
        fprintf(' example: outp = CDAA(''file_name.csv'')\n\n');
        return;
    end
    
    % Reading the file provided by its path
    [agi, t, g_raw] = ReadExprFile(tfp);
    
    % Checking whether time points values are integers
    if sum(abs(round(t) - t)) > 0
        fprintf('\nThe algorithms accepts time points as integers only!\n');
        fprintf(' Try to convert time point values into minutes\n\n');
        return;
    end
    
    % Number of genes
    n_genes = size(g_raw,1);
    % Number of time points
    n_tp = length(t);
    % Creating labels for time points to use for plotting
    lbls = cell(1,n_tp);
    for n = 1:n_tp
        lbls{n} = num2str(t(n));
    end
    
    % Asking for a file with transcription factors
    if nargin < 2
        fprintf('\nSet the path to the transcription factor agis file\n');
        fprintf(' as the second parameter\n');
        fprintf(' example: outp = CDAA(''file_name.csv'', ''tf_file name.csv'')\n');
        fprintf(' or set [] as the second parameter\n');
        fprintf(' example: outp = CDAA(''file_name.csv'',[])\n\n');
        return;
    end
    
    % Labeling Transcription Factors
    if isempty(tffp)
        tfs = true(n_genes,1);
    else
        tfs = AssignTFs(agi, tffp);
    end
    
    % Calculating normalized expression
    g = NormExpr(g_raw);
    
    % =========  Stage Separation ===========
    % Calculating normalized changes in expression
    s = NormChange(g, t);
    outp = s;
    
    if nargin < 3
        % Finding maximum values
        max_s = max(abs(s),[],2);
        % Calculating cardinalities
        G = sum(abs(s) == max_s(:,ones(1,n_tp-1)));
        % Printing the histogram
        HistPlot(G,lbls);
        % Suggesting Primary Response stage border
        %    and number of intervals
        b = find(G == max(G));
        c = 1;
        outp = [b,c];
        fprintf('\nSuggested border for the Primary Response stage\n');
        fprintf(' is b = %d and number of intervals is c = %d.\n', b, c);
        fprintf(' These or corrected values can be used as the next\n');
        fprintf(' input to the function: outp = CDAA(...,[b,c])\n\n');
        return;
    else
        % Assigning border value from the input
        b = brdr(1);
        % Assigning number of intervals value from the input
        if length(brdr) < 2
            c = 1;
        else
            c = brdr(2);
        end
    end
    
    % =========  Gene to Stage Assignment ===========
    stage_lbls = {'Initiation',...
                  'Primary Response',...
                  'Secondary Response'};
    % Initialize stage assignment
    stage = ones(n_genes,1);
    
    % Suggesting default values for cluster numbers to try
    if nargin < 4
        % in case one of the stages is missing
        if (b < 2) || (b + c >= n_tp)
            fprintf('\nSet the next input value to 0\n');
            outp = 4;
        else
            fprintf('\nSet the next input value to [0,0]\n');
            outp = [4,4];
        end
        fprintf(' to start the number of clusters selection process\n');
        fprintf(' or set the desired number of clusters\n')
        fprintf(' for stage(s) instead of zero(es)\n\n');
        return;
    end
    
    % Cluster number selection procedure
    if nargin < 5
        % Setting a default value if 0 is provided
        if n_clst(1) < 1
            n_clst(1) = 4;
        end
        % First interval is in Response stage
        if b < 2
            n_tpts = b + c;
            outp = cell(1,1);
        % Last interval is in Response stage
        elseif b + c >= n_tp
            n_tpts = b;
            outp = cell(1,1);
        else
            n_tpts = b;
            outp = cell(1,2);
        end  
        % Clustering
        [idx, cntrs] = StageCluster(g, n_tpts, n_clst(1));
        % Cluster plots
        ClusterPlot(g, n_tpts, idx, tfs, lbls);
        outp{1} = cntrs;
        fprintf('\nSet the next input value to\n');
        fprintf(' to the output of this function\n')
        fprintf(' when you are satisfied with the number of clusters\n\n');
        return;
    else
        outp = cl_cntrs;
        % Searching for the least active cluster in the first stage
        %   Clustering the first stage based on provided centroids
        n_tpts = size(cl_cntrs{1},2);
        idx = StageCluster(g, n_tpts, [], cl_cntrs{1});
        % Searching for the cluster with smallest avg diviation from 0
        sumd = inf;
        cl_n = 0;
        for n = 1:max(idx)
            sumn = norm(CenterExpr(g(idx == n, 1:n_tpts)),Inf);
            if sumn < sumd
                sumd = sumn;
                cl_n = n;
            end
        end
        % Assigning cluster with the smallest activity to the next stage
        stage = stage + (idx == cl_n);
        n_tpts = b + c;
         
        % For all 3 stages present
        if length(cl_cntrs) > 1
            % Selecting clusters for the Primary response stage
            if isempty(cl_cntrs{2})        
                % Setting a default value if none provided
                if n_clst(2) < 1
                    n_clst(2) = 4;
                end
                % Clustering
                [idx, cntrs] = StageCluster(g(stage==2,:),...
                               n_tpts, n_clst(2));
                % Cluster plots
                ClusterPlot(g(stage==2,:), n_tpts, idx, tfs(stage==2), lbls, cl_n);
                outp{2} = cntrs;
                fprintf('\nSet the next input value to\n');
                fprintf(' to the output of this function\n')
                fprintf(' when you are satisfied with the number of clusters\n\n');
            % Subclustering the Primary response stage
            else
                % Searching for the least active cluster in the second stage
                % Clustering the stage based on provided centroids
                idx = StageCluster(g(stage==2,:), n_tpts, [], cl_cntrs{2});
                idx_all = zeros(n_genes,1);
                idx_all(stage==2) = idx;
                % Searching for the cluster with smallest avg diviation from 0
                sumd = inf;
                cl_n = 0;
                for n = 1:max(idx_all)
                    sumn = norm(CenterExpr(g(idx_all == n, 1:n_tpts)),Inf);
                    if sumn < sumd
                        sumd = sumn;
                        cl_n = n;
                    end
                end
                stage = stage + (idx_all == cl_n);
            end
        end
    end
    
    % Asking for an agi of a gene of interest
    if nargin < 6
        fprintf('\nStage separation is now completed.\n');
        fprintf(' Set the agi of the gene of interest as next parameter\n');
        fprintf(' in format: CDAA(..., ''AtXgXXXXX'')\n\n');
        return;
    else
        % Finding line number that correspond to the selected gene
        goi_id = find(strcmpi(agi{1},goi));
        fprintf('\nGene with agi %s is', goi);
        outp = goi_id;
        % If no matches for agi found
        if isempty(goi_id)
            fprintf(' not found in the supplied dataset!\n\n');
            return;
        end 
        
        if ~tfs(goi_id)
            fprintf(' not');
        end
        fprintf(' assumed to be a Transcription Factor.\n');
        fprintf(' %s was classified to the %s stage.\n\n', goi,...
                stage_lbls{stage(goi_id) + (b < 2)});
    end
    
    
    % =========  Interaction inference ===========    
    % Normalizing changes in expression
    sn = s./(max(abs(s),[],2)*ones(1,n_tp - 1));
    
    % Asking for the type of regulatory connection for the gene
    if nargin < 7
        fprintf('Select whether you are looking\n');
        fprintf(' for ''regulators'' or ''targets'' of %s\n', goi);
        fprintf(' by specifying next variable.\n');
        fprintf(' For example: CDAA(..., ''regulators'')\n\n');
        outp = 'regulators';
        return;
    else
        % Specifying thresholds if not set
        if nargin < 8
            thrs = [0.4, 0.2, 0.4];
        end
        % Specifying parameters for interaction inference with regulators
        if lower(rt(1)) == 'r'
            if stage(goi_id) == 1
                fprintf('Not possible to find regulators\n');
                fprintf('  since the gene is in the earliest stage\n\n');
                return;
            elseif stage(goi_id) == 2 
                % Identifying candidate regulators
                candidates = (stage == 1) & tfs;
                % Limiting time course to 2 stages
                sn = sn(:,1:b+c-1);
                t = t(1:b+c);
                % Number of intervals
                nint = b-1;            
            else 
                % Identifying candidate regulators
                candidates = (stage == 2) & tfs;
                % Limiting time course to 2 stages
                sn = sn(:,b:n_tp-1);
                t = t(b:n_tp);
                % Number of intervals
                nint = c; 
            end
            rt_lbl = 'Regulators';
        % Specifying parameters for interaction inference with targets
        elseif lower(rt(1)) == 't'
            if ~tfs(goi_id)
                fprintf('Not possible to find targets\n');
                fprintf(' since the gene is not assumed\n');
                fprintf(' to be a Transcription Factor\n\n');
                return;
            end
            if stage(goi_id) == max(stage)
                fprintf('Not possible to find targets\n');
                fprintf(' since the gene is in the latest stage\n\n');
                return;
            elseif stage(goi_id) == 1 
                % Identifying candidate targets
                candidates = (stage == 2);
                % Limiting time course to 2 stages
                sn = sn(:,1:b+c-1);
                t = t(1:b+c);
                % Number of intervals
                nint = b-1;            
            else 
                % Identifying candidate targets
                candidates = (stage == 3);
                % Limiting time course to 2 stages
                sn = sn(:,b:n_tp-1);
                t = t(b:n_tp);
                % Number of intervals
                nint = c; 
            end
            rt_lbl = 'Targets';
        else
            fprintf('Unknown type of interaction\n\n');
            return;
        end

        % Building a dissimilarity table
        [disst, int_type, dT, cand] = GetDT(sn, t, nint,...
                                       goi_id, candidates, thrs(1), rt);
        % Table plot
        fig_n = 2;
        figure(fig_n)
        ttl = ['Dissimilarities for ' goi ' at no thr.'];
        PlotDT(disst, int_type, dT, agi{1}(cand), ttl);

        % Constructing thresholded versions of expression patterns
        if length(thrs) > 1
            % Creating a pool of putative regulators/targets and type of
            %   interaction (0 - activation, 1 - inhibition)
            reg_put = [find(cand) int_type];
            % Looking at the thresholded versions
            for thr = thrs(2:end)
                % Thresholding
                sn_thr = (-1)*(sn < -thr) + (sn > thr);
                % Building a dissimilarity table
                [disst, int_type, dT, cand] = GetDT(sn_thr,...
                                               t, nint, goi_id,...
                                               candidates, thrs(1), rt);
                % Table plot
                fig_n = fig_n + 1;
                figure(fig_n)
                ttl = ['Dissimilarities for ' goi ' at thr. = ' num2str(thr)];
                PlotDT(disst, int_type, dT, agi{1}(cand), ttl);
                % Adding new regulators/targets to the pool
                reg_put = [reg_put; [find(cand) int_type]];
            end
            % Looking for regulators/targets satisfying majority vote
            %   principle. 
            %  The number to satisfy 
            pts = length(thrs)/2;
            % The search
            reguls = [];
            for rg = unique(reg_put(:,1)')
                % Types of interaction for a connection
                int_type_rg = unique(reg_put(reg_put(:,1) == rg,2));
                % If number of thresholded patterns is not reached
                %   or more than one type of interaction found, 
                %   the connection is discarded
                if (sum(reg_put(:,1) == rg) > pts) && (length(int_type_rg) == 1)
                    reguls = [reguls; [rg int_type_rg]];
                end
            end
        else
            reguls = [find(cand) int_type];
        end

        % Printing the resulting predictions
        fprintf('%s predicted for %s:\n', rt_lbl, goi);
        outp = cell(size(reguls,1),2);
        int_type_lbl = {'activator', 'inhibitor'};
        for n = 1:size(reguls,1)
            outp{n,1} = agi{1}{reguls(n,1)};
            outp{n,2} = int_type_lbl{2 - reguls(n,2)};
            fprintf(' %s - %s\n', outp{n,1},outp{n,2});
        end
        fprintf('\n');
    end
end

function [agi, t, expr] = ReadExprFile(file_path)
%ReadExprFile extracts data from transcriptomic data file
%   agi - a cell array of gene identifiers (first column)
%   t - an array of time points (first row)
%   expr - a matrix of expression values
%   file_path - path to the .csv file with transcriptome data

    % Opening file for reading
    fileID = fopen(file_path,'r');
    % Reading the first line with sampling time points
    head = fgetl(fileID);
    head = textscan(head,'%s','delimiter',',');
    head = head{1};
    % Number of time points
    ntp = numel(head) - 1;
    
    % Reading other lines that start with genes' agis
    formatStr = ['%s', repmat('%f', 1, ntp)];
    dt = textscan(fileID,formatStr, 'Delimiter', ',',...
                  'EndOfLine', '\r\n');
    fclose(fileID);
              
    % Extracting agis
    agi = dt(:,1);
    % Initializing time points array
    t = zeros(1,ntp);
    % Initializing expression values matrix
    expr = zeros(size(dt{1},1),ntp);
    % Sorting columns in case sampling points are given out of order
    for n = 2:(ntp+1)
        t(n-1) = str2double(head(n));
        expr(:,n-1) = dt{:,n};
    end
    sorted = sortrows([t;expr]',1)';
    % Assigning time points array values
    t = sorted(1,:);
    % Assigning gene expression matrix values
    expr = sorted(2:end,:);
end

function tfs = AssignTFs(agi, tf_file_path)
%AssignTFs assigns genes with transcription factor (TF) function 
%   tfs         - trnscription factors identifier column
%   agi         - a cell array of gene identifiers
%   tf_file_path - path to the .csv file with a list of TFs 

    % Reading the file line by line
    fileID = fopen(tf_file_path,'r');
    tf_lst = textscan(fileID,'%s', 'EndOfLine', '\r\n');
    fclose(fileID);

    % Initializing the column of TF flags
    tfs = false(size(agi, 1),1);
    % Assigning 1's to the lines corresponding to TFs
    for n = 1:size(tf_lst{1},1)
        tfs = tfs | strcmpi(agi{1},tf_lst{1}(n));
    end
end

function g = CenterExpr(g_raw)
%CenterExpr removes the mean value from expression patterns 
%   g       - matrix of expression patterns with removed mean
%   g_raw   - matrix of expression patterns
    g = g_raw - mean(g_raw,2)*ones(1,size(g_raw,2));
end

function g = NormExpr(g_raw)
%NormExpr normalizes expression patterns by centering and scaling
%   g       - matrix of normalized expression patterns
%   g_raw   - matrix of expression patterns

    % Centering
    g = CenterExpr(g_raw);
    % Scaling by a standard deviation
    g = g./(std(g_raw,0,2)*ones(1,size(g,2)));
end

function s = NormChange(g, t)
%NormChange normalizes differential expression patterns by scaling
%   s - matrix of normalized differences in expression
%   g - matrix of normalized expression patterns
%   t - array of time points

    % Differences in expression
    dg = diff(g,1,2);
    % Time point differences
    dt = ones(size(g,1),1)*diff(t);
    % Getting scaled differences
    s = dg./dt;
end

function HistPlot(G,lbls)
%HistPlot plots a histogram of cardinalities
%   G       - array of cardinalities
%   lbls    - x-axis labels

    % Creating a new figure
    figure(1);
    clf;
    % Number of bars
    lg = length(G);
    % Histogram plot
    bar(1:lg,G,'BarWidth',1);
    % Labeling bars
    for r = 1:6
        text(r-0.15, G(r)+0.01*sum(G), ['$\mathcal{G}_' int2str(r) '$'],...
            'Interpreter', 'Latex');
    end
    % Setting labels
    xlim(0.5 + [0 lg])
    set(gca, 'XTick', 0.5:(lg + 0.5));
    set(gca, 'XTickLabel', lbls);
    grid on
    title('Highest expression change');
    xlabel('t');
    ylabel('Cardinality');
end

function [idx, cntrs] = StageCluster(g, n_tpts, k, c)
%StageCluster clusters genes at a certain stage
%   idx     - a column of gene cluster numbers
%   cntrs   - a matrix of resulting cluster centroids
%   g       - matrix of expression patterns
%   n_tpts  - number of time points to limit expression patterns
%   k       - number of clusters
%   c       - cluster centroids (if known)

    % Centering expression patterns
    gs = CenterExpr(g(:,1:n_tpts));
    if nargin < 4
        % Making a 1000 runs to choose the tightest clusters result
        osumd = +inf; % overall sum o squared distances
        for n = 1:1000
            [id, c, sumd] = kmeans(gs,k);
            if sum(sumd) < osumd
                idx = id;
                cntrs = c;
            end
        end
    else
        % Clustering based on known centroids
        [idx] = kmeans(gs, [], 'Start', c);
        cntrs = c;
    end
end

function ClusterPlot(g, n_tpts, idx, tfs, lbls, h_cl_n)
%ClusterPlot plots each cluster expression patterns in a separate window 
%   g       - matrix of expression patterns
%   n_tpts  - number of time points to limit expression patterns
%   idx     - a column of gene cluster numbers
%   tfs     - a coulms of flags of whether the gene is a transcription factor
%   lbls    - x-axis labels
%   h_cl_n  - parent cluster number (in case of subclustering)

    % Taking the mean out
    gs = CenterExpr(g(:,1:n_tpts));
    fg = 1;
    % Going through each cluster
    for n = 1:length(unique(idx))
        fg = fg + 1;
        figure(fg);
        clf;
        % Selecting genes that belong to the current cluster
        gsn = gs(idx == n, :);
        % Cluster plot
        plot(1:n_tpts, gsn, 'Color', 0.65*ones(1,3));
        % Ploting transcription factor patterns with different color
        if sum((idx == n) & tfs) > 0
            hold on
            gsn_tf = gs((idx == n) & tfs, :);
            plot(1:n_tpts, gsn_tf, 'b');
        end
        grid on;
        % Setting labels
        set(gca, 'XTick', 1:length(lbls));
        set(gca, 'XTickLabel', lbls);
        xlabel('t');
        if nargin < 6
            cl_n = num2str(n);
        else
            cl_n = [num2str(h_cl_n) '.' num2str(n)];
        end
        title(['Cluster ' cl_n ...
               ' (n_{genes} = ' num2str(sum(idx == n))...
               ', n_{TF} = ' num2str(sum((idx == n) & tfs)) ')']);
    end
end

function [disst, int_type, dT, cnd] = GetDT(sn, t, nint, goi_id,...
                                               candidates, thr, rt)
%GetDT calculates dissimilarity table
%   disst       - dissimilarity table
%   int_type    - types of interaction (inhibition/activation)
%   dT          - greatest common divisior for time intervals
%   cnd         - updated candidates set
%   sn          - normalized changes in expression
%   t           - set of time points for regulator's and target's stages
%   nint        - number of intervals in regulator's stage
%   goi_id      - row of sn that corresponds to the gene
%   candidates  - genes classified to the preceeding stage
%   thr         - threshold for minimal dissimilarity
%   rt          - relationship type (regulator/target)

    % Finding the greatest common divisior for time intervals
    %   Getting time intervals
    intls = diff(t);
    %   The search
    dT = min(intls);
    for intl = intls
        if gcd(dT,intl) < dT
            dT = gcd(dT,intl);
        end
    end
    
    % Zeroth-order approximation
    %   Number of copies for each time interval
    nc = intls/dT;
    %   The approximation
    sn0 = [];
    for n = 1:size(sn,2)
        sn0 = [sn0, repmat(sn(:,n), 1, nc(n))];
    end
    
    % Number of time intervals in the regulator's stage
    nti_reg = sum(nc(1:nint));
    % Number of time intervals in the target's stage
    nti_tgt = sum(nc) - nti_reg;
    
    % Building a dissimilarity table
    %   Number of time shifts
    n_sh = min(nti_reg,nti_tgt);
    %   Window size
    wnd = sum(nc) - n_sh;
    %   Dissimilarities for activation and inhibition effects
    disa = []; disi = [];
    
    % When looking for regulators
    if lower(rt(1)) == 'r'
        %   Candidate ragulators' patterns portion
        cand_prt = sn0(candidates,1:wnd);
        %   Target's pattern copies
        targ_ptrn = sn0(goi_id*ones(1,sum(candidates)),:);

        %   Alignment
        for n = 0:min(nti_reg,nti_tgt)
            % Target's portion
            targ_prt = targ_ptrn(:,n+1:n+wnd);
            % Dissimilarities
            disa = [disa, mean(abs(cand_prt - targ_prt),2)];
            disi = [disi, mean(abs(cand_prt + targ_prt),2)];
        end
    % When looking for targets
    else
        %   Ragulator's pattern copies
        reg_prt = sn0(goi_id*ones(1,sum(candidates)),1:wnd);
        %   Candidate target's patterns
        cand_ptrn = sn0(candidates,:);

        %   Alignment
        for n = 0:min(nti_reg,nti_tgt)
            % Target's portion
            cand_prt = cand_ptrn(:,n+1:n+wnd);
            % Dissimilarities
            disa = [disa, mean(abs(reg_prt - cand_prt),2)];
            disi = [disi, mean(abs(reg_prt + cand_prt),2)];
        end
    end
    %   Interaction types
    int_types = disa < disi;
    %   Dissimilarity table
    disst = disa.*int_types + disi.*(1-int_types);
    %   Number of dissimilarities
    n_diss = size(disst,2);
    % Deleting entrances where minimum dissimilarity
    %   is at no shift or outside Resulator's stage region
    % Minimum dissimilarity
    min_diss = min(disst,[],2);
    % Candidates to delete
    dlt_lns = ~((min_diss == min(disst(:,2:max(2,end-1)),[],2)) &...
                (min_diss < thr));
    % Removing candidates
    cnd = candidates;
    cnd(candidates) = ~dlt_lns;
    disst(dlt_lns,:) = [];
    int_types(dlt_lns,:) = [];
    min_diss(dlt_lns) = [];
    if n_diss > 2
        disst(:,[1,n_diss]) = [];
        int_types(:,[1,n_diss]) = [];
        n_diss = n_diss - 2;
        % Choosing an interaction type
        int_type = int_types(:,1);
        % Searching for an interaction type at minimum dissimilarity
        for n = 2:n_diss
            % Minimum dissimilarity lines
            mdlns = (disst(:,n) == min_diss);
            % correcting the interaction type
            int_type(mdlns) = int_types(mdlns, n);
        end
    else
        disst(:,1) = [];
        int_type = int_types(:,2);
    end
end

function PlotDT(disst, int_type, dT, agis, ttl)
%PlotDT plots dissimilarity table
%   disst       - dissimilarity table
%   int_type    - types of interaction (inhibition/activation)
%   dT          - greatest common divisior for time intervals
%   agis        - set of genes' agis
%   ttl         - title (optional)

    % Clearing the figure
    clf;
    % Interaction type labels
    itt = ['a', 'i'];
    % Number of regulators
    n_regs = size(disst,1);
    % Range of time shifts
    xrange = dT:dT:dT*size(disst,2);
    % Table plot
    imagesc(xrange, [], disst);
    
    % Setting labels
    set(gca,'XTick', xrange);
    set(gca,'CLim', [0 1]);
    set(gca,'YTick', 1:size(disst,1))
    yt_lbl = cell(1,n_regs);
    for n = 1:n_regs
        yt_lbl{n} = [agis{n} ' - ' itt(2 - int_type(n))];
    end
    set(gca,'YTickLabel',yt_lbl)
    xlabel('Delay');
    if nargin > 4
        title(ttl);
    end
    colorbar;
end
