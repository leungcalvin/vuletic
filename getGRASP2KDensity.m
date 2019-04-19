function [grid,totalChargeDensity] = getGRASP2KDensity(Npoints,orbitalList,occupationNumbers);
Norbs = numel(occupationNumbers);

totalChargeDensity = zeros(Npoints,1);
for idx=1:Norbs
    disp(orbitalList(idx))
    grasptable = table2array(import_rwf(string(orbitalList(idx)), 1, Npoints));
    r = grasptable(:,1);large  = grasptable(:,2);small  = grasptable(:,3);
    totalChargeDensity = totalChargeDensity+ occupationNumbers(idx) .* (large.^2 + small.^2);
    grid = r;
end
end