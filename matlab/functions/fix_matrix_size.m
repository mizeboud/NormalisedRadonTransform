function fixed_matrix = fix_matrix_size(change_mat,ref_mat)
% -- update size of im_resz to actual data crevSig
while any(size(change_mat) ~= size(ref_mat))  
    if size(change_mat,1) <  size(ref_mat,1)
        change_mat(end:size(ref_mat,1),:) = nan;
    elseif size(change_mat,2) < size(ref_mat,2)
        change_mat(:,end:size(ref_mat,2)) = nan;
    elseif size(change_mat,1) > size(ref_mat,1)
        change_mat = change_mat(1:size(ref_mat,1),:) ;
    elseif size(change_mat,2) > size(ref_mat,2)
        change_mat = change_mat(:,1:size(ref_mat,2)) ;
    end
end
fixed_matrix = change_mat;
end
    