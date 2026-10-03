%% sim_Stance_Swing_Biped.m
%
% Sagittal-plane BIPED gait: two point-foot legs, each the same 2-link leg
% used in sim_Leg_Swing.m, hanging from a common pelvis that carries the
% HAT point mass (mH). Leg 2 runs the same cycle as leg 1, half a cycle later.
%
% Method (same idea as the single-leg scripts: prescribe the motion, then
% inverse dynamics):
%   1. Prescribe the PELVIS path and the two FOOT paths (task space).
%   2. Inverse kinematics gives each leg's hip angle and knee angle.
%   3. Each leg's inverse dynamics (swing model + moving-hip term + foot
%      force) gives the joint torques.
%   4. Total ground force from whole-body momentum, shared between the two
%      feet by a smooth load-transfer weight (double support).
%
% Angle convention (same as sim_Leg_Swing.m animation):
%   th  = absolute thigh angle from vertical, forward (+x) positive
%   kap = knee flexion (relative), shank angle = th - kap
%   u(a) = [sin a; -cos a] is the unit vector at angle a (0 = straight down)
% -------------------------------------------------------------------

clear; clc; close all;

%% 1. Parameters ---------------------------------------------------------
P.mH  = 44.0;                         % kg   HAT (point mass at pelvis)
P.m1  = 7.0;   P.m2  = 4.5;           % kg   thigh, shank+foot
P.L1  = 0.45;  P.L2  = 0.45;          % m
P.lc1 = 0.20;  P.lc2 = 0.19;          % m    COM distances
P.I1  = 0.10;  P.I2  = 0.05;          % kg m^2
P.g   = 9.81;

P.Ts = 0.6;                           % s  stance duration of one leg
P.Tw = 0.4;                           % s  swing duration of one leg
P.S    = 0.40;                        % m  step length (footprint to opposite footprint)
P.ybar = 0.80;                        % m  mean pelvis height
P.A    = 0.03;                        % m  peak-to-peak vertical pelvis oscillation
P.clr  = 0.06;                        % m  swing-foot clearance

P.T     = P.Ts + P.Tw;                % stride period
P.Tstep = P.T/2;                      % one step
P.Td    = P.Ts - P.Tstep;             % double-support time per step
P.v     = P.S/P.Tstep;                % forward speed
assert(P.Td > 0 && P.Ts > 2*P.Td, 'Need 0 < Td and Ts > 2*Td');

nCycles = 4;
dt = 0.005;
t  = 0:dt:nCycles*P.T;
N  = numel(t);

fprintf('Stride T = %.2f s, step %.2f s, double support %.2f s per step, speed %.2f m/s\n', ...
    P.T, P.Tstep, P.Td, P.v);

selfTest(P);

%% 2. Prescribed motion -> IK -> inverse dynamics ------------------------
pel   = zeros(2,N);
th    = zeros(2,N);   kap = zeros(2,N);      % (leg, time)
tau   = zeros(2,2,N);                        % (joint hip/knee, leg, time)
lam   = zeros(2,2,N);                        % (x/y, leg, time) ground force on each foot
lamT  = zeros(2,N);                          % total ground force
wgt   = zeros(2,N);                          % load share of each leg
stance= false(2,N);
kneeW = zeros(2,2,N);  footW = zeros(2,2,N); % (x/y, leg, time) world positions
ext   = zeros(2,N);                          % leg extension d/(L1+L2)

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
    % total ground force from whole-body momentum (like groundReactionForce)
    lamT(:,k) = P.mH*(pk.pdd + [0; P.g]) + F;
    for i = 1:2
        lam(:,i,k) = wgt(i,k)*lamT(:,k);                     % load sharing
        tau(:,i,k) = tau0(:,i) - Jf(:,:,i)' * lam(:,i,k);   % joint torques
    end
end

tauHip  = squeeze(tau(1,:,:));   % 2 x N  (row = leg)
tauKnee = squeeze(tau(2,:,:));
trunkRes = -(tauHip(1,:) + tauHip(2,:));   % torque the trunk would have to supply

%% 3. Checks -------------------------------------------------------------
mg = (P.mH + 2*(P.m1+P.m2))*P.g;
fprintf('\n--- Checks ---\n');
fprintf('Mean vertical ground force: %.1f N   (body weight %.1f N)\n', mean(lamT(2,:)), mg);
fprintf('Max leg extension: %.1f %% of leg length (keep < ~97 %%)\n', 100*max(ext(:)));
fprintf('Min vertical force per foot: %.2f N (should be >= 0)\n', min(min(squeeze(lam(2,:,:)))));
mu = 0;
for i = 1:2
    ly = squeeze(lam(2,i,:))';  lx = squeeze(lam(1,i,:))';
    m = ly > 50;  mu = max([mu, abs(lx(m)./ly(m))]);
end
fprintf('Max |Fx/Fy| on a loaded foot: %.2f (friction coefficient needed)\n', mu);
fprintf('Max |tau_hip| = %.0f N.m, max |tau_knee| = %.0f N.m\n', max(abs(tauHip(:))), max(abs(tauKnee(:))));
fprintf('Max |trunk residual| = %.0f N.m (hip torques should sum to 0 for a free trunk)\n', max(abs(trunkRes)));
footY = squeeze(footW(2,:,:));            % 2 x N
fprintf('Lowest foot height: %.4f m, stance-foot height error: %.1e m\n', ...
    min(footY(:)), max(abs(footY(stance))));

%% 4. Plots --------------------------------------------------------------
cLeg = [0 0.447 0.741; 0.850 0.325 0.098];

figure('Name','Biped: joint angles and torques','Color','w');
subplot(2,1,1); hold on; grid on;
plot(t, rad2deg(th(1,:)),  '-',  'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, rad2deg(kap(1,:)), '--', 'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, rad2deg(th(2,:)),  '-',  'Color',cLeg(2,:), 'LineWidth',1.5);
plot(t, rad2deg(kap(2,:)), '--', 'Color',cLeg(2,:), 'LineWidth',1.5);
xlabel('time (s)'); ylabel('angle (deg)');
legend('hip 1','knee 1','hip 2','knee 2','Location','best');
title('Joint angles (leg 1 blue, leg 2 red)');

subplot(2,1,2); hold on; grid on;
plot(t, tauHip(1,:),  '-',  'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, tauKnee(1,:), '--', 'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, tauHip(2,:),  '-',  'Color',cLeg(2,:), 'LineWidth',1.5);
plot(t, tauKnee(2,:), '--', 'Color',cLeg(2,:), 'LineWidth',1.5);
xlabel('time (s)'); ylabel('torque (N\cdot m)');
legend('\tau_{hip} 1','\tau_{knee} 1','\tau_{hip} 2','\tau_{knee} 2','Location','best');
title('Required joint torques');

figure('Name','Biped: ground reaction forces','Color','w');
subplot(2,1,1); hold on; grid on;
plot(t, squeeze(lam(2,1,:)), 'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, squeeze(lam(2,2,:)), 'Color',cLeg(2,:), 'LineWidth',1.5);
plot(t, lamT(2,:), 'k', 'LineWidth',1.2);
yline(mg,'k--','body weight');
xlabel('time (s)'); ylabel('vertical force (N)');
legend('foot 1','foot 2','total','Location','best');
title('Vertical ground reaction force');

subplot(2,1,2); hold on; grid on;
plot(t, squeeze(lam(1,1,:)), 'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, squeeze(lam(1,2,:)), 'Color',cLeg(2,:), 'LineWidth',1.5);
xlabel('time (s)'); ylabel('horizontal force (N)');
legend('foot 1','foot 2','Location','best');
title('Horizontal ground reaction force');

figure('Name','Biped: support and trunk residual','Color','w');
subplot(2,1,1); hold on; grid on;
plot(t, wgt(1,:), 'Color',cLeg(1,:), 'LineWidth',1.5);
plot(t, wgt(2,:), 'Color',cLeg(2,:), 'LineWidth',1.5);
xlabel('time (s)'); ylabel('load share'); ylim([-0.05 1.05]);
legend('leg 1','leg 2','Location','best');
title('Load transfer (both > 0 = double support)');
subplot(2,1,2);
plot(t, trunkRes, 'k', 'LineWidth',1.5); grid on;
xlabel('time (s)'); ylabel('N\cdot m');
title('Trunk residual  -(\tau_{hip1}+\tau_{hip2})  (diagnostic)');

%% 5. Animation ----------------------------------------------------------
animateBiped(t, pel, kneeW, footW, lam, P, 5);


% =========================================================================
%                          LOCAL FUNCTIONS
% =========================================================================

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
% With a fixed hip and no foot force, the new leg dynamics must reproduce the
% original swing-model torque M*qdd + C*qd + G.
    worst = 0;
    for n = 1:20
        q = [0.8*randn; 0.1 + rand]; qd = randn(2,1); qdd = 3*randn(2,1);
        t1 = leg_inverse_dynamics(q, qd, qdd, [0;0], P);
        [M,C,G] = twoLinkDynamicsMatrices_swing(q(1), q(2), qd(1), qd(2), P);
        t2 = M*qdd + C*qd + G;
        worst = max(worst, max(abs(t1 - t2)));
    end
    fprintf('Self-test (fixed hip vs. swing model): max torque difference %.1e N.m\n', worst);
end

function animateBiped(t, pel, kneeW, footW, lam, P, slow)
    fps = 60;  dtt = t(2) - t(1);
    cLeg = [0 0.447 0.741; 0.850 0.325 0.098];
    fmax = max(hypot(lam(1,:,:), lam(2,:,:)), [], 'all');
    fscale = 0.4/fmax;

    fig = figure('Name','Biped Gait Animation','Color','w');
    ax = axes(fig); hold(ax,'on'); axis(ax,'equal'); ylim(ax,[-0.1 1.5]);
    xGround = [-5, pel(1,end)+5];
    plot(ax, xGround, [0 0], 'k-', 'LineWidth', 2);
    tx = ceil(xGround(1)/0.5)*0.5 : 0.5 : xGround(2);
    plot(ax, tx, zeros(size(tx)), 'k|', 'MarkerSize', 10);
    xlabel(ax,'x (m, walking direction)'); ylabel(ax,'height (m)');

    trunk = plot(ax,[0 0],[0 0],'-','Color',[0.4 0.4 0.4],'LineWidth',6);
    legL  = gobjects(2,1);  arr = gobjects(2,1);
    for i = 1:2
        legL(i) = plot(ax,[0 0 0],[0 0 0],'-o','Color',cLeg(i,:),'LineWidth',3, ...
            'MarkerFaceColor',cLeg(i,:),'MarkerSize',6);
        arr(i) = quiver(ax,0,0,0,0,'AutoScale','off','Color',cLeg(i,:), ...
            'LineWidth',2,'MaxHeadSize',0.8);
    end
    pelPt  = plot(ax,0,0,'ks','MarkerFaceColor','k','MarkerSize',9);
    titleH = title(ax,'');

    N = numel(t);
    clk = tic;
    while true
        if ~isvalid(fig), return; end
        el = toc(clk)/slow;
        done = el >= t(end);
        if done, k = N; else, k = floor(el/dtt) + 1; end

        pp = pel(:,k);
        set(trunk,'XData',[pp(1) pp(1)],'YData',[pp(2) pp(2)+0.5]);
        set(pelPt,'XData',pp(1),'YData',pp(2));
        for i = 1:2
            kn = kneeW(:,i,k);  ft = footW(:,i,k);  F = lam(:,i,k);
            set(legL(i),'XData',[pp(1) kn(1) ft(1)],'YData',[pp(2) kn(2) ft(2)]);
            set(arr(i),'XData',ft(1),'YData',ft(2),'UData',F(1)*fscale,'VData',F(2)*fscale);
        end
        xlim(ax, pp(1) + [-1 1]*1.2);
        set(titleH,'String',sprintf('t = %.2f s   |R_1| = %.0f N   |R_2| = %.0f N', ...
            t(k), hypot(lam(1,1,k),lam(2,1,k)), hypot(lam(1,2,k),lam(2,2,k))));
        drawnow limitrate;
        if done, break; end
    end
end