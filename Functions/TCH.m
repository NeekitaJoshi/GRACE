%
% Here we try to solve the system (5) for R (covariance matrix) as in
% Ferreira et al. (2016), so one should used the variances in the diagonal,
% that is, r11, r22, r33, r44. For the standard deviation, don't forget the
% sqrt(r11), etc.

% Reference:

% Ferreira VG, Montecino HDC, Yakubu CI, Heck B (2016) Uncertainties of the
% Gravity Recovery and Climate Experiment time-variable gravity-field
% solutions based on three-cornered hat method. J Appl Remote Sens 10:015015.
% doi: 10.1117/1.JRS.10.015015

function [S,R] = TCH(x1,x2,x3,x4)
% The "reference" can be anyone of the four inputs, here we choose the
% last one (i.e. x4). See Eq. (2) of Ferreira et al. (2016).
y14 = x1-x4; y24 = x2-x4; y34 = x3-x4;

% The estimates of measurement (co-)variances denoted by sij as in Eq. (2)
% of Ferreira et al. (2016). Thanks to Matlab we have:
S = cov([y14 y24 y34]);
% for input the co-variances in Eqs. (8) and (9), we need them as a vector
% for convenience.
s = [S(1,1) S(1,2) S(1,3) S(2,2) S(2,3) S(3,3)];
% 3-vector need in Eqs. (5)-(7).
u = [1; 1; 1];
% Initial conditions as suggested at the relations (10)
r14 = 0; r24 = 0; r34 = 0;
r44 = (2*u'*inv(S)*u)^-1;
x0 = [r14 r24 r34 r44];
% Here we use the matlab optmization toolbox for computing the N-free
% parameters.
opts = optimset('Algorithm','active-set','TolX',1e-10,'TolCon',1e-10,'Display','off');
x = fmincon(@(r) myfun(r,s),x0,[],[],[],[],[],[],@(r) mycon(r,s),opts);
% Once the free-parameters have been determined, we can compute the remain
% elements of R using as in Eq. (7), which provides
r14 = x(1); r24 = x(2); r34 = x(3); r44 = x(4);
r11 = s(1)-r44+2*r14;
r12 = s(2)-r44+r14+r24;
r13 = s(3)-r44+r14+r34;
r22 = s(4)-r44+2*r24;
r23 = s(5)-r44+r24+r34;
r33 = s(6)-r44+2*r34;
% Finaly the full Allan's covariance matrix is
R = [r11 r12 r13 r14; r12 r22 r23 r24; r13 r23 r33 r34; r14 r24 r34 r44];
end