function MI_mat_perm = create_MI_mat_perm(sig_trim, sim_num)

MI_mat_perm = nan([size(sig_trim, 2), size(sig_trim, 2), sim_num]);
for i = 1 : sim_num
    if mod(i, 100) == 0
        disp(i)
    end
    sig_trim_perm = sig_trim;
    perm_dist = randperm(size(sig_trim, 1), size(sig_trim_perm, 2));    
    for j = 1 : size(sig_trim_perm, 2)
        sig_trim_perm(:, j) = ...
            circshift(sig_trim_perm(:, j), perm_dist(j));
    end
    
    MI_mat_perm(:, :, i) = calculate_MI_mat(sig_trim_perm);
end
end