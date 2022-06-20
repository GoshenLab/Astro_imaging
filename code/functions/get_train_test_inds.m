function [train_ind, test_ind] = ...
    get_train_test_inds(train_prop, total_length)


train_length = round(train_prop * total_length);
train_start = randi(total_length);
train_ind =...
    train_start : ...
    mod((train_start + train_length), total_length) - 1;
if isempty(train_ind)
    
    if (mod((train_start + train_length), total_length) - 1) < 0
        % if the test is exactly the first frames
        
        train_ind =...
            train_start : total_length; 
        % rounding issues fix:
        if length(train_ind) > train_length
            train_ind =...
                (train_start+1) : total_length;
        elseif length(train_ind) < train_length
             train_ind =...
                (train_start-1) : total_length;
        end
    else
        train_ind = ...
            [train_start : total_length, ...
            1 : mod((train_start + train_length), total_length) - 1];
    end
end

test_vec = ones(total_length, 1);
test_vec(train_ind) = 0;
test_ind = find(test_vec);

end
