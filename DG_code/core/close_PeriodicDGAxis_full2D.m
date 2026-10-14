function ops = close_PeriodicDGAxis_full2D(ops)
%CLOSE_PERIODICDGAXIS_FULL2D Treat the last/first DG faces as one interface.
% The central and dissipative parts use the same flux as interior faces.
n = ops.nDof;
left = sparse(1,1,1,1,n);
right = sparse(1,n,1,1,n);
ops.H1 = ops.H1 - 0.5*ops.leftLift*right + 0.5*ops.rightLift*left;
seam = 0.5*(ops.leftLift-ops.rightLift)*(left-right);
ops.H2Interior = ops.H2Interior + seam;
ops.H2Boundary = sparse(n,n);
ops.H2 = ops.H2Interior;
end
