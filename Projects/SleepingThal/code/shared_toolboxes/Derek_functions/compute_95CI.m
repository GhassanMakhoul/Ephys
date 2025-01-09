function [mean_data,CI,lower_CI,higher_CI,N] = compute_95CI(input_data,dim)
% Compute the 95% confidence interval of data
%   Can take in an input of any shape. Use the dim parameter to indicate
%   which dimension to take the 95
% Inputs:
%   input_data -> a n-dimensional matrix to take the 95CI of. 
%   dim -> dimension to take the 95CI over
% Outputs:
%   mean_data -> means
%   CI -> CIs
%   lower_CI -> lower CI
%   higher_CI -> higher CI
%   N -> num samples in each

if sum(size(input_data)==1) == (length(size(input_data))-1)
    % size(input_data,1) == 1 | size(input_data,2) == 1
    cleaned_data = input_data(~isnan(input_data));
    mean_data = mean(cleaned_data,'all','omitnan');
    
    CI = tinv([0.975],length(cleaned_data)-1) .* (std(cleaned_data,0,'all','omitnan')) ./ sqrt(length(cleaned_data));
    
    lower_CI = mean_data - CI;
    higher_CI = mean_data + CI;
else
    num_dims = ndims(input_data);
    curr_dims = [1:num_dims];
    curr_dims(dim) = [];    
    rearranged_input_data = permute(input_data,[dim ,curr_dims]);
    

    N = sum(~isnan(rearranged_input_data),1);

    mean_data = squeeze(mean(rearranged_input_data,1,'omitnan'));
    SEM = (std(rearranged_input_data,0,1,'omitnan')) ./ sqrt(N);
    fixed_N_minus_1 = N - 1;
    fixed_N_minus_1(fixed_N_minus_1<0) = 0;
    CI_standard = tinv([0.975],fixed_N_minus_1);
    CI = squeeze(CI_standard .* SEM);

    lower_CI = mean_data - CI;
    higher_CI = mean_data + CI;

end