% Here is the objective function as per Eq. (8) in Ferreira et al. (2016)
function F = myfun(r,s)
F = (r(1)^2 + r(2)^2 + r(3)^2 + (r(1) + r(2) - r(4) + s(2))^2 + (r(1) +...
    r(3) - r(4) + s(3))^2 + (r(2) + r(3) - r(4) + s(5))^2)/(- s(6)*s(2)^2 +...
    2*s(2)*s(3)*s(5) - s(4)*s(3)^2 - s(1)*s(5)^2 + s(1)*s(4)*s(6))^(2/3);
end