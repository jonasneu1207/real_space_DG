function result = collect_sweep_qtbm2D(points,mat,options,status)
%COLLECT_SWEEP_QTBM2D Common output for fixed and self-consistent sweeps.
count = numel(points); current = nan(1,count); density = current; biases = current;
valid = false(1,count);
for k = 1:count
    p = points{k}; biases(k) = p.bias; valid(k) = p.converged;
    if valid(k)
        current(k) = p.transport.j; density(k) = p.transport.currentDensity;
    end
end
result = struct('status',status,'converged',all(valid),'options',options, ...
    'mat',mat,'points',{points},'sweep',struct('bias',biases,'current',current, ...
    'currentDensity',density,'converged',valid,'requestedBias',options.bias));
end
