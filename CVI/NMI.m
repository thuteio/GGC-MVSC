function nmi_score = NMI(true_labels, pred_labels)
% 快速版本的NMI，适用于数值型标签
    
    if ~isnumeric(true_labels) || ~isnumeric(pred_labels)
        error('此快速版本仅适用于数值标签');
    end
    
    n = length(true_labels);
    
    % 使用accumarray快速计算联合分布
    [true_unique, ~, true_idx] = unique(true_labels);
    [pred_unique, ~, pred_idx] = unique(pred_labels);
    
    contingency = accumarray([true_idx, pred_idx], 1);
    
    % 计算概率分布
    joint_prob = contingency / n;
    true_prob = sum(joint_prob, 2);
    pred_prob = sum(joint_prob, 1);
    
    % 计算互信息（向量化计算）
    [i_idx, j_idx] = find(joint_prob > 0);
    p_ij = joint_prob(joint_prob > 0);
    p_i = true_prob(i_idx);
    p_j = pred_prob(j_idx)';
    
    mi = sum(p_ij .* log(p_ij ./ (p_i .* p_j)));
    
    % 计算熵
    h_true = -sum(true_prob(true_prob > 0) .* log(true_prob(true_prob > 0)));
    h_pred = -sum(pred_prob(pred_prob > 0) .* log(pred_prob(pred_prob > 0)));
    
    % 计算NMI
    if h_true == 0 && h_pred == 0
        nmi_score = 1;
    elseif h_true == 0 || h_pred == 0
        nmi_score = 0;
    else
        nmi_score = mi / ((h_true + h_pred) / 2);
    end
end