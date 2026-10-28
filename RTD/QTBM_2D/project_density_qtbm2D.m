function [nPoisson,info] = project_density_qtbm2D(transport,mat)
%PROJECT_DENSITY_QTBM2D Conservative piecewise-linear nodal dual-cell average.
% End X cells have half width; final Y column is a periodic duplicate.
n = transport.n;
if ~isequal(size(n),[mat.Nx,mat.Ny]) || any(~isfinite(n(:))) || any(n(:)<0) ...
        || max(abs(n-n(:,1)),[],'all')>1e-10*max(1,max(n(:)))
    error('QTBM2D:Density','Expected nonnegative, finite, Y-homogeneous point density.');
end
x = mat.x(:)*1e-9; dx = diff(x);
left = [0;dx]/2; right = [dx;0]/2; volume = left+right;
line = n(:,1);
integral = left.*(3*line+[line(1);line(1:end-1)])/4 ...
    +right.*(3*line+[line(2:end);line(end)])/4;
nPoisson = repmat(integral./volume,1,mat.Ny);
width = (mat.y(end)-mat.y(1))*1e-9;
inputCharge = trapz(x,line)*width; outputCharge = sum(integral)*width;
info = struct('inputParticlesPerLength',inputCharge, ...
    'outputParticlesPerLength',outputCharge, ...
    'relativeConservationError',abs(outputCharge-inputCharge)/max(abs(inputCharge),realmin), ...
    'xControlVolumeM',volume,'reconstruction','piecewise-linear point density');
end
