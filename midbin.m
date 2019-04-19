function dw = midbin(r)
dw_left = [0; diff(r)];
dw_right = [diff(r); 0];
dw = (dw_left + dw_right) * 0.5;
    % zero on the left of the list
    end