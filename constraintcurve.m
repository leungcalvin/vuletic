function curves = constraintcurve(psX,psY,FYFX,sensitivity_Hz)
    curves = sensitivity_Hz./(abs(psY .* ones(size(FYFX)) - psX .* FYFX));
    %[curve] = NxM; N different particle masses and M different values of
    %FYFX
end