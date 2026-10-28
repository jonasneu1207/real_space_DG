function [nodes,weights] = gauss_rule_qtbm2D(order)
%GAUSS_RULE_QTBM2D Gauss-Legendre nodes and weights on [-1,1].
persistent orders nodeCache weightCache
index = find(orders==order,1);
if isempty(index)
    k = (1:order-1).'; beta = k./sqrt(4*k.^2-1);
    [V,D] = eig(diag(beta,1)+diag(beta,-1));
    [nodes,indexSort] = sort(diag(D));
    weights = 2*V(1,indexSort).'.^2;
    orders(end+1) = order; index = numel(orders);
    nodeCache{index} = nodes; weightCache{index} = weights;
else
    nodes = nodeCache{index}; weights = weightCache{index};
end
end
