function ramp = bias_ramp_qtbm2D(mat)
%BIAS_RAMP_QTBM2D Same initial active-region drop as the DG RTD.
active = mat.x(end)-mat.x(1)-mat.Lc1-mat.Lc2;
if active<=0, error('QTBM2D:Geometry','No active interval between contacts.'); end
ramp = min(1,max(0,(mat.x(:)-mat.x(1)-mat.Lc1)/active));
end
