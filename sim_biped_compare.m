%% sim_Biped_Compare.m   (Q3, part c)
%
% Compare the joint torques and ground reaction forces of three men with body
% masses 70, 85 and 100 kg, walking with the biped model of sim_Biped_Gait.m
% (part b). The simulation itself (prescribed pelvis/foot paths -> inverse
% kinematics -> inverse dynamics -> ground force) is the part-b code, wrapped
% in runBiped() so it can be run once per person. Nothing in it was changed.
%
% PARAMETER SOURCE (one source for every body parameter)
%   Winter, D. A. (2009). Biomechanics and Motor Control of Human Movement,
%   4th ed. Hoboken, NJ: Wiley.  doi:10.1002/9780470549148
%   Openly accessible copy of the pages used (Fig. 4.1 and Table 4.1, book pp. 83-86):
%   https://courses.grainger.illinois.edu/me481/sp2021/Anthro-Winter.pdf
%   (University of Illinois, ME 481 course reading; accessed 3 Oct 2026)
%   Table 4.1 (adult-male data of Dempster 1955), M = body mass, L = segment length:
%       thigh                  mass 0.100 M   COM 0.433 L from hip    r_g(COM) 0.323 L
%       "foot and leg"         mass 0.061 M   COM 0.606 L from knee   r_g(COM) 0.416 L
%       HAT (head, arms, trunk) mass 0.678 M  (point mass at the pelvis, as in part b)
%   Fig. 4.1 (Drillis & Contini 1966), H = stature:
%       hip height 0.530 H, knee height 0.285 H, ankle height 0.039 H
%       => thigh L1 = 0.245 H, knee-to-ground L2 = 0.285 H (point foot),
%          knee-to-ankle length of the "foot and leg" segment = 0.246 H.
%   Moment of inertia about each segment's own COM:  I = m (r_g L)^2.
%
% ASSUMPTIONS
%   * The problem fixes only the body masses, so all three men get the same
%     stature H (edit the vector H below to give each man his own).
%   * The gait pattern is that of part (b), with every length scaled by
%     leg length / 0.90 m (geometric similarity); step and stride TIMES are unchanged.
%   * With equal H the kinematics are identical, so every torque and force
%     scales exactly with M. The normalised plots therefore collapse onto one
%     curve (a built-in check); with different H they separate slightly.
% -------------------------------------------------------------------

clear; clc; close all;

M = [70 85 100];                 % kg
H = [1.75 1.75 1.75];            % m  stature (assumed, see above)
nCycles = 2;
dt = 0.005;

nS = numel(M);
for s = 1:nS
    P(s) = buildParams(M(s), H(s));
    fprintf('--- %d kg, H = %.2f m: L1 = %.3f m, L2 = %.3f m, mH = %.1f kg, m1 = %.1f kg, m2 = %.2f kg\n', ...
        M(s), H(s), P(s).L1, P(s).L2, P(s).mH, P(s).m1, P(s).m2);
    selfTest(P(s));
    R(s) = runBiped(P(s), nCycles, dt);
end

%% Metrics ---------------------------------------------------------------
Mt = struct('FyPk',[], 'FxPk',[], 'hipPk',[], 'kneePk',[], 'mu',[], 'ext',[], 'FyMin',[]);
for s = 1:nS
    r  = R(s);
    Fy = squeeze(r.lam(2,:,:));   Fx = squeeze(r.lam(1,:,:));      % 2 x N (leg, time)
    Mt.FyPk(s)   = max(Fy(:));
    Mt.FyMin(s)  = min(Fy(:));
    Mt.FxPk(s)   = max(abs(Fx(:)));
    Mt.hipPk(s)  = max(abs(r.tauHip(:)));
    Mt.kneePk(s) = max(abs(r.tauKnee(:)));
    Mt.ext(s)    = 100*max(r.ext(:));
    mu = 0;
    for i = 1:2
        m = Fy(i,:) > 50;
        mu = max([mu, abs(Fx(i,m)./Fy(i,m))]);
    end
    Mt.mu(s) = mu;
end
BW = M*9.81;

row = @(label, vals, fmt) fprintf(['%-40s' repmat(fmt, 1, nS) '\n'], label, vals);
fprintf('\n%-40s', 'Peak values over the walk');
for s = 1:nS, fprintf('%12s', sprintf('%d kg', M(s))); end
fprintf('\n');
row('Peak vertical ground force (N)',        Mt.FyPk,           '%12.0f');
row('Peak vertical ground force (body wt)',  Mt.FyPk./BW,       '%12.2f');
row('Peak |horizontal| ground force (N)',    Mt.FxPk,           '%12.0f');
row('Peak |hip torque| (N.m)',               Mt.hipPk,          '%12.0f');
row('Peak |knee torque| (N.m)',              Mt.kneePk,         '%12.0f');
row('Peak |hip torque| per kg (N.m/kg)',     Mt.hipPk./M,       '%12.2f');
row('Peak |knee torque| per kg (N.m/kg)',    Mt.kneePk./M,      '%12.2f');
row('Ratio to first man: mass',              M/M(1),            '%12.3f');
row('Ratio to first man: peak vertical force', Mt.FyPk/Mt.FyPk(1), '%12.3f');
row('Ratio to first man: peak hip torque',   Mt.hipPk/Mt.hipPk(1),   '%12.3f');
row('Ratio to first man: peak knee torque',  Mt.kneePk/Mt.kneePk(1), '%12.3f');
fprintf('\nSanity checks: min vertical force on a foot >= 0 -> %s | max leg extension %.1f %% | friction needed (|Fx/Fy|) %.2f\n', ...
    mat2str(all(Mt.FyMin >= -1e-6)), max(Mt.ext), max(Mt.mu));

%% Plots -----------------------------------------------------------------
compareFigure('Biped: torques and ground forces (absolute)',   R, P, M, false);
compareFigure('Biped: torques and ground forces (per kg / per body weight)', R, P, M, true);

figure('Name','Biped: peak values vs body mass','Color','w');
subplot(1,3,1); plot(M, Mt.FyPk,   'o-', 'LineWidth',1.5); grid on;
xlabel('body mass (kg)'); ylabel('N');     title('Peak vertical ground force');
subplot(1,3,2); plot(M, Mt.hipPk,  'o-', 'LineWidth',1.5); grid on;
xlabel('body mass (kg)'); ylabel('N\cdot m'); title('Peak |hip torque|');
subplot(1,3,3); plot(M, Mt.kneePk, 'o-', 'LineWidth',1.5); grid on;
xlabel('body mass (kg)'); ylabel('N\cdot m'); title('Peak |knee torque|');


% =========================================================================
%                          LOCAL FUNCTIONS
% =========================================================================

function P = buildParams(M, H)
% Body parameters of a man of mass M (kg) and stature H (m) from Winter (2009),
% plus the part-(b) gait pattern scaled with leg length.
    P.mH  = 0.678*M;                  % HAT
    P.m1  = 0.100*M;                  % thigh
    P.m2  = 0.061*M;                  % foot and leg (shank + foot)
    P.L1  = 0.245*H;                  % hip -> knee
    P.L2  = 0.285*H;                  % knee -> ground (point foot)
    Lka   = 0.246*H;                  % knee -> ankle (length used by Table 4.1)
    P.lc1 = 0.433*P.L1;               % thigh COM from hip
    P.lc2 = 0.606*Lka;                % foot-and-leg COM from knee
    P.I1  = P.m1*(0.323*P.L1)^2;      % about own COM
    P.I2  = P.m2*(0.416*Lka)^2;
    P.g   = 9.81;

    k = (P.L1 + P.L2)/0.90;           % leg length relative to the part-(b) leg (0.90 m)
    P.Ts = 0.6;  P.Tw = 0.4;          % same timing as part (b)
    P.S    = 0.40*k;                  % step length
    P.ybar = 0.80*k;                  % mean pelvis height
    P.A    = 0.03*k;                  % vertical pelvis oscillation
    P.clr  = 0.06*k;                  % swing-foot clearance

    P.T     = P.Ts + P.Tw;
    P.Tstep = P.T/2;
    P.Td    = P.Ts - P.Tstep;
    P.v     = P.S/P.Tstep;
    assert(P.Td > 0 && P.Ts > 2*P.Td, 'Need 0 < Td and Ts > 2*Td');
end

function R = runBiped(P, nCycles, dt)
% The simulation loop of sim_Biped_Compare.m (section 2), unchanged, returned as a struct.
    t  = 0:dt:nCycles*P.T;
    N  = numel(t);
    pel   = zeros(2,N);
    th    = zeros(2,N);   kap = zeros(2,N);
    tau   = zeros(2,2,N);
    lam   = zeros(2,2,N);
    lamT  = zeros(2,N);
    wgt   = zeros(2,N);
    stance= false(2,N);
    kneeW = zeros(2,2,N);  footW = zeros(2,2,N);
    ext   = zeros(2,N);

    for k = 1:N
        [pk, leg] = gaitTargets(t(k), P);
        pel(:,k) = pk.p;
        F = zeros(2,1);  tau0 = zeros(2,2);  Jf = zeros(2,2,2);
        for i = 1:2
            [q, qd, qdd, ~, d] = leg_IK(leg(i).f - pk.p, leg(i).fd - pk.pd, leg(i).fdd - pk.pdd, P);
            [tau0(:,i), Fl, Jf(:,:,i)] = leg_inverse_dynamics(q, qd, qdd, pk.pdd, P);
            F = F + Fl;
            th(i,k) = q(1);  kap(i,k) = q(2);  ext(i,k) = d/(P.L1+P.L2);
            wgt(i,k) = leg(i).w;  stance(i,k) = leg(i).stance;
            [kn, ft] = legFK(q(1), q(2), P);
            kneeW(:,i,k) = pk.p + kn(:);   footW(:,i,k) = pk.p + ft(:);
        end
        lamT(:,k) = P.mH*(pk.pdd + [0; P.g]) + F;
        for i = 1:2
            lam(:,i,k) = wgt(i,k)*lamT(:,k);
            tau(:,i,k) = tau0(:,i) - Jf(:,:,i)' * lam(:,i,k);
        end
    end

    R.t = t;  R.pel = pel;  R.th = th;  R.kap = kap;
    R.tauHip  = squeeze(tau(1,:,:));      % 2 x N  (row = leg)
    R.tauKnee = squeeze(tau(2,:,:));
    R.lam = lam;  R.lamT = lamT;  R.wgt = wgt;  R.stance = stance;
    R.kneeW = kneeW;  R.footW = footW;  R.ext = ext;
end

function compareFigure(figName, R, P, M, normalise)
% Leg 1 over one stride (leg 2 is the same curve half a stride later).
% 0-60 % of the stride is stance, 60-100 % is swing.
    c = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
    mask = R(1).t <= P(1).T + 1e-9;
    x = 100*R(1).t(mask)/P(1).T;
    names = arrayfun(@(m) sprintf('%d kg', m), M, 'UniformOutput', false);
    lw   = [6 3.5 1.5];            % normalised plot: thick, medium, thin so overlapping
    lsty = {'-', '--', '-'};       % curves stay visible (they coincide when all heights are equal)
    ttl = {'Hip torque (leg 1)', 'Knee torque (leg 1)', ...
           'Vertical ground force (foot 1)', 'Horizontal ground force (foot 1)'};
    if normalise
        yl = {'N\cdot m / kg', 'N\cdot m / kg', 'body weights', 'body weights'};
    else
        yl = {'N\cdot m', 'N\cdot m', 'N', 'N'};
    end

    figure('Name', figName, 'Color', 'w');
    for q = 1:4
        subplot(2,2,q); hold on; grid on;
        for s = 1:numel(R)
            switch q
                case 1, y = R(s).tauHip(1,:);
                case 2, y = R(s).tauKnee(1,:);
                case 3, y = squeeze(R(s).lam(2,1,:))';
                case 4, y = squeeze(R(s).lam(1,1,:))';
            end
            if normalise
                if q <= 2, y = y/M(s); else, y = y/(M(s)*P(s).g); end
            end
            if normalise
                plot(x, y(mask), lsty{s}, 'Color', c(s,:), 'LineWidth', lw(s));
            else
                plot(x, y(mask), 'Color', c(s,:), 'LineWidth', 1.5);
            end
        end
        xline(100*P(1).Ts/P(1).T, 'k:');
        xlabel('% of stride'); ylabel(yl{q}); title(ttl{q});
        if q == 1, legend(names, 'Location', 'best'); end
    end
end

% ---------------- part (b) functions, unchanged ---------------------------

function [s,sd,sdd] = minjerk(t,T)
    tau = min(max(t,0),T)/T;
    s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;
    sd  = (30*tau.^2 - 60*tau.^3 + 30*tau.^4)/T;
    sdd = (60*tau - 180*tau.^2 + 120*tau.^3)/T^2;
end

function s = smoothstep5(u)
% 0 -> 1 with zero slope and zero curvature at both ends
    u = min(max(u,0),1);
    s = 6*u^5 - 15*u^4 + 10*u^3;
end

function [pel, leg] = gaitTargets(t, P)
% Pelvis and foot paths (world frame) with velocities and accelerations.
% Leg i has its own cycle clock; leg 2 is shifted by T/2.
    w  = 2*pi/P.Tstep;
    ph = w*(t - P.Td/2);                  % lowest point in the middle of double support
    pel.p   = [P.v*t;  P.ybar - (P.A/2)*cos(ph)];
    pel.pd  = [P.v;    (P.A/2)*w*sin(ph)];
    pel.pdd = [0;      (P.A/2)*w^2*cos(ph)];

    off = [0, P.T/2];
    for i = 1:2
        s   = t + off(i);
        tau = s - floor(s/P.T)*P.T;       % time within this leg's cycle, 0 = heel-strike
        tHS = t - tau;                    % heel-strike time of the current cycle
        Fa  = P.v*(tHS + P.Ts/2);         % footprint: pelvis is over it at mid-stance

        if tau < P.Ts                     % ---- stance: foot planted
            leg(i).stance = true;
            leg(i).f = [Fa; 0];  leg(i).fd = [0; 0];  leg(i).fdd = [0; 0];
        else                              % ---- swing: footprint -> footprint + 2S
            leg(i).stance = false;
            ts = tau - P.Ts;
            [s1,s1d,s1dd] = minjerk(ts, P.Tw);
            a = pi*ts/P.Tw;  uu = sin(a);  cc = cos(a);
            x = Fa + 2*P.S*s1;   xd = 2*P.S*s1d;   xdd = 2*P.S*s1dd;
            y   = P.clr*uu^4;                                   % zero slope and curvature at both ends
            yd  = P.clr*4*pi*uu^3*cc/P.Tw;
            ydd = P.clr*pi^2*(12*uu^2*cc^2 - 4*uu^4)/P.Tw^2;
            leg(i).f = [x; y];  leg(i).fd = [xd; yd];  leg(i).fdd = [xdd; ydd];
        end

        % load share: ramp up after heel-strike, ramp down before toe-off
        if tau < P.Td
            leg(i).w = smoothstep5(tau/P.Td);
        elseif tau < P.Ts - P.Td
            leg(i).w = 1;
        elseif tau < P.Ts
            leg(i).w = 1 - smoothstep5((tau - (P.Ts-P.Td))/P.Td);
        else
            leg(i).w = 0;
        end
    end
end

function [q, qd, qdd, J, d] = leg_IK(r, rd, rdd, P)
% Foot position relative to the hip -> [th; kap], with rates.
% r = L1*u(th) + L2*u(th-kap)
    d = norm(r);
    if d >= P.L1 + P.L2
        error('Leg fully stretched (d = %.3f). Lower ybar or step length S.', d);
    end
    kap   = acos(max(-1, min(1, (d^2 - P.L1^2 - P.L2^2)/(2*P.L1*P.L2))));
    alpha = atan2(r(1), -r(2));
    del   = atan2(P.L2*sin(kap), P.L1 + P.L2*cos(kap));
    th    = alpha + del;
    ths   = th - kap;

    u1 = [sin(th);  -cos(th)];   e1 = [cos(th);  sin(th)];
    u2 = [sin(ths); -cos(ths)];  e2 = [cos(ths); sin(ths)];
    J  = [P.L1*e1 + P.L2*e2,  -P.L2*e2];       % dr/d[th; kap]

    qd   = J \ rd(:);
    thsd = qd(1) - qd(2);
    qdd  = J \ (rdd(:) + P.L1*u1*qd(1)^2 + P.L2*u2*thsd^2);
    q    = [th; kap];
end

function [tau0, Fleg, J] = leg_inverse_dynamics(q, qd, qdd, ap, P)
% Joint torques [tau_hip; tau_knee] of one leg whose hip accelerates with ap,
% WITHOUT any foot force (tau0). The foot force enters afterwards as
%     tau = tau0 - J' * lambda.
% Fleg = sum of m*(a_COM + g) over the leg's bodies (force the hip + foot
% must supply). With ap = 0 this is exactly M*qdd + C*qd + G of sim_Leg_Swing.m.
    th = q(1);  ths = q(1) - q(2);
    thsd = qd(1) - qd(2);   thsdd = qdd(1) - qdd(2);
    u1 = [sin(th);  -cos(th)];   e1 = [cos(th);  sin(th)];
    u2 = [sin(ths); -cos(ths)];  e2 = [cos(ths); sin(ths)];

    rho1dd = P.lc1*e1*qdd(1) - P.lc1*u1*qd(1)^2;                      % COM accel relative to hip
    rho2dd = P.L1*e1*qdd(1) - P.L1*u1*qd(1)^2 ...
           + P.lc2*e2*thsdd - P.lc2*u2*thsd^2;
    gv = [0; P.g];
    F1 = P.m1*(ap(:) + rho1dd + gv);
    F2 = P.m2*(ap(:) + rho2dd + gv);

    Jc1 = [P.lc1*e1, [0;0]];
    Jc2 = [P.L1*e1 + P.lc2*e2,  -P.lc2*e2];
    tau0 = Jc1'*F1 + Jc2'*F2 + [P.I1*qdd(1) + P.I2*thsdd; -P.I2*thsdd];

    Fleg = F1 + F2;
    J = [P.L1*e1 + P.L2*e2,  -P.L2*e2];                               % foot Jacobian
end

function [knee, foot] = legFK(th, kap, P)
% Foot and knee relative to the hip (hip at origin, x forward, y up).
    d1 = [sin(th), -cos(th)];
    ck = cos(kap); sk = sin(kap);
    d2 = [ck*d1(1)+sk*d1(2), -sk*d1(1)+ck*d1(2)];
    knee = P.L1*d1;
    foot = knee + P.L2*d2;
end

function [M, C, G] = twoLinkDynamicsMatrices_swing(theta_h, theta_k, thetadot_h, thetadot_k, p)
% Unchanged from sim_Leg_Swing.m (used only in the self-test).
    m1=p.m1; m2=p.m2; L1=p.L1; lc1=p.lc1; lc2=p.lc2; I1=p.I1; I2=p.I2; g=p.g;
    M11 = I1+I2+L1^2*m2+2*L1*lc2*m2*cos(theta_k)+lc1^2*m1+lc2^2*m2;
    M12 = -I2-L1*lc2*m2*cos(theta_k)-lc2^2*m2;
    M22 = I2+lc2^2*m2;
    M = [M11 M12; M12 M22];
    h = L1*lc2*m2*sin(theta_k);
    C = [-h*thetadot_k,  h*(thetadot_k-thetadot_h);
          h*thetadot_h,  0];
    G1 = g*((L1*m2+lc1*m1)*sin(theta_h) + lc2*m2*sin(theta_h-theta_k));
    G2 = -g*lc2*m2*sin(theta_h-theta_k);
    G = [G1; G2];
end

function selfTest(P)
% With a fixed hip and no foot force, the leg dynamics must reproduce the
% original swing-model torque M*qdd + C*qd + G.
    worst = 0;
    for n = 1:20
        q = [0.8*randn; 0.1 + rand]; qd = randn(2,1); qdd = 3*randn(2,1);
        t1 = leg_inverse_dynamics(q, qd, qdd, [0;0], P);
        [M,C,G] = twoLinkDynamicsMatrices_swing(q(1), q(2), qd(1), qd(2), P);
        t2 = M*qdd + C*qd + G;
        worst = max(worst, max(abs(t1 - t2)));
    end
    fprintf('    self-test (fixed hip vs. swing model): max torque difference %.1e N.m\n', worst);
end