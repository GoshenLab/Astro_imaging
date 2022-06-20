function MI_mat = calculate_MI_mat(sig)

sig_length = size(sig, 1);
sig_num = size(sig, 2);

MI_mat = nan(sig_num);

for i = 1 : sig_num
    for j = i + 1 : sig_num
        
        x = sig(:, i);
        y = sig(:, j);
        
        J = [sum(~x & ~y), sum(~x & y); ...
            sum(x & ~y), sum(x & y)] ./ sig_length;
        
        % this is to fix for 0 prob events. The diag yields the self info.
        J(J==0) = nan;
        
        % mutual information in bits:
        
        MI = sum(sum(J .* log2(J ./ ...
            (sum(J, 2, 'omitnan') * sum(J, 1, 'omitnan'))), 'omitnan'),...
            'omitnan');
        
        
        MI_mat(i, j) = MI;
        MI_mat(j, i) = MI;
    end
    
end
end