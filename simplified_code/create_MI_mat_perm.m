function MI_mat_perm = create_MI_mat_perm(sig, sim_num)

MI_mat_perm = nan([size(sig, 2), size(sig, 2), sim_num]);
for i = 1 : sim_num
    if mod(i, 100) == 0
        disp(i)
    end
    sig_perm = sig;
    perm_dist = randperm(size(sig, 1), size(sig_perm, 2));    
    for j = 1 : size(sig_perm, 2)
        sig_perm(:, j) = ...
            circshift(sig_perm(:, j), perm_dist(j));
    end
    
    MI_mat_perm(:, :, i) = calculate_MI_mat(sig_perm);
end
end