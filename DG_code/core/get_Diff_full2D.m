function [A, rhs, info] = get_Diff_full2D(p)
%GET_DIFF_FULL2D Transport/flux scaffold for full physical 2D Wigner-DG.
%
% This function does not assemble the global 4D transport matrix yet. It
% exposes the rectangular tensor-product building blocks that will later
% enter the X- and Y-transport operators by Kronecker products.

A = [];
rhs = [];
info = struct;
info.fullMatrixAssembled = false;
info.reason = 'Boundary/flux scaffold only; no global 4D matrix is built here.';
info.dofOrder = p.index.order;
info.localMassXYSize = size(p.dg.MXY);
info.localDXSize = size(p.dg.DX);
info.localDYSize = size(p.dg.DY);
info.relativeAxSize = size(p.relative.Ax);
info.relativeAySize = size(p.relative.Ay);
info.faces.leftTraceSize = size(p.dg.faces.local.left.trace);
info.faces.rightTraceSize = size(p.dg.faces.local.right.trace);
info.faces.bottomTraceSize = size(p.dg.faces.local.bottom.trace);
info.faces.topTraceSize = size(p.dg.faces.local.top.trace);
end
