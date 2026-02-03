function sys = normalizeParameters(sys, value)
%NORMALIZEPARAMETERS 
% Torben Warnecke - 11/06/2024

HphiEq = sys.H.phi.equality.c + sys.H.phi.equality.t - sys.H.phi.equality.f;
HphiIneq = sys.H.phi.inequality.c + sys.H.phi.inequality.t - sys.H.phi.inequality.f;

for k = 1:size(HphiEq,1)
    a = sum(HphiEq(k,:).^2);
    HphiEq(k,:) = HphiEq(k,:) .* sqrt(value/a);
    %HphiEq(k,:) = HphiEq(k,:) / max(abs(HphiEq(k,:))) * value;
end
for k = 1:size(HphiIneq,1)
    a = sum(HphiIneq(k,:).^2);
    HphiIneq(k,:) = HphiIneq(k,:) .* sqrt(value/a);
    %HphiIneq(k,:) = HphiIneq(k,:) / max(abs(HphiIneq(k,:))) * value;
end

sys.H.phi.equality = HphiEq;
sys.H.phi.inequality = HphiIneq;
end

