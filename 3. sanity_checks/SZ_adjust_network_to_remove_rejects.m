%% Aim: for the rejected signal, replace amplitude and latencies of N1 and N2 with NaN. 
% Would be used during sanity checks (correlation between amplitude and latencies & correlation between amplitude and distance)
function out = SZ_adjust_network_to_remove_rejects(out)

% get num of electrodes for this patient 
num_elecs = size(out.elecs, 2);

% Pre-allocate n1_adj and n2_adj matrix with nan
for row = 1:num_elecs
    if size(out.elecs(row).avg, 1) >= 1
        out.elecs(row).n1_adj = zeros(num_elecs, 2); 
        out.elecs(row).n2_adj = zeros(num_elecs, 2); 
    end
end


% for each patient's out file, loop over all electrodes 
for ich = 1:num_elecs
    
    %% adj n1
    % if the current elecs is not a stimulating elecs, skip it as it
    % contains no data to be altered.
    if isempty(out.elecs(ich).arts), continue; end
        
    % extract n1 
    n1 = out.elecs(ich).N1(:,1:2);
    
    %{
    % sz commented out on 2025/04/03
    % replace the 0s with NaN in the new_version of N1 amp and latencies
    for jch = 1:num_elecs
        if n1(jch, 1) == 0 && n1(jch, 2) == 0 % if both amp and lat are zeros
            n1(jch, 1) = NaN;
            n1(jch, 2) = NaN;
        end
    end
    %}
   

    % exclude the n1 if the ccep signal are rejected during pipeline
    % get the keep for the current stimulating electrode
    keep_vector = out.rejection_details(1).reject.keep(ich,:);
    % identify indices where keep_vector is NaN or 0
    invalid_indices = isnan(keep_vector) | (keep_vector == 0);
    % replace corresponding elements in n1 amplitude and latency with zeros
    n1(invalid_indices, 1) = 0;
    n1(invalid_indices, 2) = 0;
        
    % fill n1_adj with the adjusted n1 amp and lat
    out.elecs(ich).n1_adj = n1;


    %% adj n2 
    % extract n1 
    n2 = out.elecs(ich).N2(:,1:2);
    
    %{
    % replace the 0s with NaN in the new_version of N1 amp and latencies
    for jch = 1:num_elecs
        if n2(jch, 1) == 0 && n2(jch, 2) == 0 % if both amp and lat are zeros
            n2(jch, 1) = NaN;
            n2(jch, 2) = NaN;
        end
    end
    %}

    % replace corresponding elements in n1 amplitude and latency with NaN
    n2(invalid_indices, 1) = 0;
    n2(invalid_indices, 2) = 0;
        
    % fill n1_adj with the adjusted n1 amp and lat
    out.elecs(ich).n2_adj = n2;

end

%{

%% Aim: for the rejected signal, replace amplitude and latencies of N1 and N2 with NaN. 
% Would be used during sanity checks (correlation between amplitude and latencies & correlation between amplitude and distance)
function out = SZ_adjust_network_to_remove_rejects(out)


% get num of electrodes for this patient 
num_elecs = size(out.elecs, 2);


% Pre-allocate n1_adj and n2_adj matrix with nan
for row = 1:num_elecs
    if size(out.elecs(row).avg, 1) >= 1
        out.elecs(row).n1_adj = nan(num_elecs, 2); 
        out.elecs(row).n2_adj = nan(num_elecs, 2); 
    end
end


% for each patient's out file, loop over all electrodes 
for ich = 1:num_elecs
    
    %% adj n1
    % if the current elecs is not a stimulating elecs, skip it as it
    % contains no data to be altered.
    if isempty(out.elecs(ich).arts), continue; end
        
    % extract n1 
    n1 = out.elecs(ich).N1(:,1:2);
        
    % replace the 0s with NaN in the new_version of N1 amp and latencies
    for jch = 1:num_elecs
        if n1(jch, 1) == 0 && n1(jch, 2) == 0 % if both amp and lat are zeros
            n1(jch, 1) = NaN;
            n1(jch, 2) = NaN;
        end
    end
   
    % exclude the n1 if the ccep signal are rejected during pipeline
    % get the keep for the current stimulating electrode
    keep_vector = out.rejection_details(1).reject.keep(ich,:);
    % identify indices where keep_vector is NaN or 0
    invalid_indices = isnan(keep_vector) | (keep_vector == 0);
    % replace corresponding elements in n1 amplitude and latency with NaN
    n1(invalid_indices, 1) = NaN;
    n1(invalid_indices, 2) = NaN;
        
    % fill n1_adj with the adjusted n1 amp and lat
    out.elecs(ich).n1_adj = n1;


    %% adj n2 
    % extract n1 
    n2 = out.elecs(ich).N2(:,1:2);
        
    % replace the 0s with NaN in the new_version of N1 amp and latencies
    for jch = 1:num_elecs
        if n2(jch, 1) == 0 && n2(jch, 2) == 0 % if both amp and lat are zeros
            n2(jch, 1) = NaN;
            n2(jch, 2) = NaN;
        end
    end

    % replace corresponding elements in n1 amplitude and latency with NaN
    n2(invalid_indices, 1) = NaN;
    n2(invalid_indices, 2) = NaN;
        
    % fill n1_adj with the adjusted n1 amp and lat
    out.elecs(ich).n2_adj = n2;

end


%}