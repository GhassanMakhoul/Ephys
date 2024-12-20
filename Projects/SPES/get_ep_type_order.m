function out = get_ep_type_order(in_list, orig_type, pat_orig_order)

out = nan(length(in_list),1);

for i = 1:length(in_list)
    curr_str = in_list(i);
    s = split(curr_str, "_");
    pat_id = s(1);
    idx = find(contains(pat_orig_order, pat_id));
    if length(idx) ~= 1; error("ERROR: not exactly one idx"); end
    out(i) = orig_type(idx);






end