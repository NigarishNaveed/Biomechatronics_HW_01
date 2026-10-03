%% sim_Bowling_Elbow.m
%
% Q4: Advantage a cricket bowler extracts by straightening the elbow
% during the bowling action (shoulder + elbow two-link arm).
%
% Model: planar two-link arm from the class notes (Euler-Lagrange).
%   th1 = shoulder angle, measured from the horizontal x-axis (x = direction
%         of bowling), counter-clockwise positive
%   th2 = elbow angle RELATIVE to the upper arm
%         th2 = 0 -> arm fully straight,  th2 > 0 -> elbow flexed
%         (forearm trails behind the upper arm, as in a cocked throw)
%   Equations of motion (class notes, step 8):
%       M(th) thdd + C(th,thd) thd + G(th) = [tau_shoulder; tau_elbow]
%
% Range of motion analysed (ICC):  from the horizontal crossing to the point
% of delivery.
%   horizontal crossing : arm horizontal behind the body  -> th1 = 180 deg
%   vertical crossing   : arm vertical overhead           -> th1 =  90 deg
%   delivery            : 20 deg past the vertical        -> th1 =  70 deg
% so the shoulder sweeps 180 -> 70 deg (clockwise, 110 deg).
%
% Method (same as the leg scripts): PRESCRIBE the motion, compute the
% torques by INVERSE dynamics, then feed those torques into the FORWARD
% dynamics (ode45) as a consistency check.
%
% Cases: elbow straightens by 0 (straight arm), 10, 15 and 20 deg between
% the horizontal crossing and delivery. The shoulder motion, masses and
% lengths are IDENTICAL in all cases, so any difference is due to the elbow.
% -------------------------------------------------------------------

clear; clc; close all;

%% 1. Parameters ---------------------------------------------------------
% Segment proportions follow Winter, "Biomechanics and Motor Control of
% Human Movement" (4th ed., 2009), Table 4.1 (Dempster data) and Fig. 4.1.
Mbody = 85;      % kg  bowler mass (illustrative)
Hbody = 1.85;    % m   bowler height (illustrative)

p.L1  = 0.186*Hbody;                 % upper arm length
Lf    = 0.146*Hbody;                 % forearm length
Lh    = 0.108*Hbody;                 % hand length
p.L2  = Lf + Lh;                     % elbow to fingertips (ball release point)

p.m1  = 0.028*Mbody;                 % upper arm mass
p.lc1 = 0.436*p.L1;                  % COM from shoulder
p.I1  = p.m1*(0.322*p.L1)^2;         % inertia about own COM (radius of gyration 0.322 L)

% link 2 = forearm + hand + cricket ball, combined into one rigid body
mf = 0.016*Mbody;  df = 0.430*Lf;          % forearm mass, COM from elbow
mh = 0.006*Mbody;  dh = Lf + 0.506*Lh;     % hand mass, COM from elbow
mB = 0.156;        dB = p.L2;              % ball (kg), held at the fingertips
p.m2  = mf + mh + mB;
p.lc2 = (mf*df + mh*dh + mB*dB)/p.m2;
p.I2  = mf*(0.303*Lf)^2 + mf*(df-p.lc2)^2 ...
      + mh*(0.297*Lh)^2 + mh*(dh-p.lc2)^2 ...
      + mB*(dB-p.lc2)^2;
p.g   = 9.81;

% Prescribed motion (assumptions)
m.Tb       = 0.12;    % s    time from horizontal crossing to delivery
m.th_start = 180;     % deg  horizontal crossing
m.th_rel   = 70;      % deg  delivery = 20 deg past vertical (90 - 20)
m.w0       = 6;       % rad/s shoulder speed at horizontal crossing
m.w1       = 25;      % rad/s shoulder speed at delivery
m.k        = 2.5;     % elbow extension is back-loaded: at delivery the elbow
                      % is still opening at k*(extension)/Tb  (rad/s)

cases   = [0 10 15 20];            % deg of elbow straightening (0 = straight arm)
iccLimit = 15;                     % deg, ICC permitted extension
nC      = numel(cases);
t       = linspace(0, m.Tb, 601)'; % output time grid
nT      = numel(t);

doAnimate = true;                  % animation of all four arms
slowFactor = 30;                   % animation slow-motion factor

%% 2. Simulate every case -------------------------------------------------
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
for c = 1:nC
    dl = cases(c);

    % --- prescribed motion, inverse dynamics --------------------------
    q = zeros(nT,2); qd = zeros(nT,2); qdd = zeros(nT,2); tau = zeros(nT,2);
    for k = 1:nT
        [q(k,:), qd(k,:), qdd(k,:)] = arm_desired(t(k), m, dl);
        [M,C,G] = twoLinkDynamicsMatrices_arm(q(k,:)', qd(k,:)', p);
        tau(k,:) = (M*qdd(k,:)' + C*qd(k,:)' + G)';
    end

    % --- forward simulation (consistency check) -----------------------
    tauFun = @(tt) arm_inverse_dynamics(tt, m, dl, p);
    x0 = [q(1,:)'; qd(1,:)'];
    [~, Y] = ode45(@(tt,x) twoLinkEOM_arm(tt, x, tauFun, p), t, x0, opts);
    trackErr = max(max(abs(Y(:,1:2) - q)))*180/pi;

    % --- ball (fingertip) kinematics ----------------------------------
    elbow = zeros(2,nT); tip = zeros(2,nT); vtip = zeros(2,nT);
    for k = 1:nT
        [elbow(:,k), tip(:,k), vtip(:,k)] = armKinematics(q(k,:)', qd(k,:)', p);
    end
    speed = hypot(vtip(1,:), vtip(2,:))';

    % --- results --------------------------------------------------------
    R(c).dl = dl;
    R(c).q = q; R(c).qd = qd; R(c).tau = tau;
    R(c).elbow = elbow; R(c).tip = tip; R(c).vtip = vtip; R(c).speed = speed;
    R(c).vRel   = speed(end);
    R(c).launch = atan2d(vtip(2,end), vtip(1,end));      % + above horizontal
    R(c).elbowRate = rad2deg(qd(end,2));                 % deg/s at delivery
    R(c).peakTau   = max(abs(tau));                      % [shoulder elbow]
    R(c).Wtot   = trapz(t, tau(:,1).*qd(:,1) + tau(:,2).*qd(:,2));
    R(c).Welbow = trapz(t, tau(:,2).*qd(:,2));
    R(c).trackErr = trackErr;
end

%% 3. Results table ------------------------------------------------------
v0 = R(1).vRel;
fprintf('Arm: L1 = %.3f m, L2 = %.3f m (reach %.3f m), m1 = %.2f kg, m2 = %.2f kg (incl. ball)\n', ...
    p.L1, p.L2, p.L1+p.L2, p.m1, p.m2);
fprintf('Shoulder: %d -> %d deg in %.2f s;  shoulder angle at delivery = %.1f deg\n\n', ...
    m.th_start, m.th_rel, m.Tb, rad2deg(R(1).q(end,1)));
fprintf('%-11s %-8s %-10s %-10s %-8s %-14s %-12s %-12s %-10s\n', ...
    'Extension','ICC','v_ball','gain','gain','elbow rate','peak tau_sh','peak tau_el','total work');
fprintf('%-11s %-8s %-10s %-10s %-8s %-14s %-12s %-12s %-10s\n', ...
    '(deg)','','(m/s)','(m/s)','(%)','(deg/s)','(N.m)','(N.m)','(J)');
for c = 1:nC
    if R(c).dl == 0, icc = '-'; elseif R(c).dl <= iccLimit, icc = 'legal'; else, icc = 'ILLEGAL'; end
    fprintf('%-11d %-8s %-10.2f %-10.2f %-8.1f %-14.0f %-12.0f %-12.0f %-10.1f\n', ...
        R(c).dl, icc, R(c).vRel, R(c).vRel - v0, 100*(R(c).vRel - v0)/v0, ...
        R(c).elbowRate, R(c).peakTau(1), R(c).peakTau(2), R(c).Wtot);
end
fprintf('\nBall launch angle at delivery (same for all cases): %.1f deg (negative = below horizontal)\n', R(1).launch);
fprintf('Max tracking error, forward sim vs prescribed motion: %.2e deg\n', max([R.trackErr]));
fprintf('  (tiny -> the inverse-dynamics torques reproduce the prescribed motion)\n');

i10 = find(cases==10); i15 = find(cases==15); i20 = find(cases==20);
if ~isempty(i10) && ~isempty(i15) && ~isempty(i20)
    fprintf('\n15 deg vs 10 deg: +%.2f m/s extra;   20 deg vs 15 deg: +%.2f m/s extra\n', ...
        R(i15).vRel - R(i10).vRel, R(i20).vRel - R(i15).vRel);
    fprintf('Going from the legal 15 deg to 20 deg adds %.1f %% more speed gain.\n', ...
        100*((R(i20).vRel - v0)/(R(i15).vRel - v0) - 1));
end

%% 4. Plots --------------------------------------------------------------
col = [0 0 0; 0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188];
lbl = arrayfun(@(d) sprintf('%d deg extension', d), cases, 'UniformOutput', false);
lbl{1} = 'straight arm (0 deg)';

figure('Name','Prescribed motion','Color','w');
subplot(2,1,1);
plot(t, rad2deg(R(1).q(:,1)), 'LineWidth',1.5); grid on; hold on;
yline(90,'k:','vertical'); 
xlabel('time (s)'); ylabel('\theta_1 (deg)');
title('Shoulder angle: horizontal crossing (180 deg) to delivery (70 deg) - same for all cases');
subplot(2,1,2); hold on; grid on;
for c = 1:nC
    plot(t, rad2deg(R(c).q(:,2)), 'Color',col(c,:), 'LineWidth',1.5);
end
xlabel('time (s)'); ylabel('\theta_2 (deg)');
legend(lbl,'Location','best');
title('Elbow flexion (0 = straight): arm straightens by the amount shown');

figure('Name','Ball speed','Color','w');
subplot(2,1,1); hold on; grid on;
for c = 1:nC
    plot(t, R(c).speed, 'Color',col(c,:), 'LineWidth',1.5);
end
xlabel('time (s)'); ylabel('fingertip speed (m/s)');
legend(lbl,'Location','northwest');
title('Speed of the ball (fingertip) relative to the shoulder');
subplot(2,1,2);
gain = [R.vRel] - v0;
b = bar(cases, gain); grid on; b.FaceColor = [0.2 0.5 0.8];
for c = 1:nC
    text(cases(c), gain(c), sprintf(' %.2f m/s (%.1f%%)', gain(c), 100*gain(c)/v0), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
xline(iccLimit,'r--','ICC limit 15 deg');
xlabel('elbow extension (deg)'); ylabel('release-speed gain (m/s)');
title('Advantage over a straight arm at the point of delivery');
ylim([0 1.2*max(gain)+eps]);

figure('Name','Joint torques','Color','w');
subplot(2,1,1); hold on; grid on;
for c = 1:nC, plot(t, R(c).tau(:,1), 'Color',col(c,:), 'LineWidth',1.5); end
xlabel('time (s)'); ylabel('\tau_{shoulder} (N\cdot m)'); legend(lbl,'Location','best');
title('Required shoulder torque');
subplot(2,1,2); hold on; grid on;
for c = 1:nC, plot(t, R(c).tau(:,2), 'Color',col(c,:), 'LineWidth',1.5); end
xlabel('time (s)'); ylabel('\tau_{elbow} (N\cdot m)'); legend(lbl,'Location','best');
title('Required elbow torque');

figure('Name','Arm positions','Color','w');
ax = axes; hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
idx = round(linspace(1, nT, 7));
for c = [1 i15]
    for k = idx
        plot(ax, [0 R(c).elbow(1,k) R(c).tip(1,k)], [0 R(c).elbow(2,k) R(c).tip(2,k)], ...
            '-', 'Color', col(c,:), 'LineWidth', 1.5);
    end
end
for c = 1:nC
    plot(ax, R(c).tip(1,:), R(c).tip(2,:), '-', 'Color', col(c,:), 'LineWidth', 1);
    plot(ax, R(c).tip(1,end), R(c).tip(2,end), 'o', 'Color', col(c,:), 'MarkerFaceColor', col(c,:));
end
plot(ax, 0, 0, 'ks', 'MarkerFaceColor','k');
xlabel('x (m, bowling direction)'); ylabel('y (m)');
title('Arm snapshots (black = straight, other = 15 deg) and fingertip paths');

%% 5. Animation ----------------------------------------------------------
if doAnimate
    animateArms(t, R, col, lbl, p, slowFactor);
end


% =========================================================================
%                          LOCAL FUNCTIONS
% =========================================================================

function [y,yd,ydd] = hermite5(t,T,y0,y1,v0,v1)
% Quintic from y0 to y1 over [0,T] with end velocities v0, v1 and zero
% end accelerations.
    u = min(max(t,0),T)/T;
    y   = y0*(1-10*u^3+15*u^4-6*u^5) + y1*(10*u^3-15*u^4+6*u^5) ...
        + T*v0*(u-6*u^3+8*u^4-3*u^5) + T*v1*(-4*u^3+7*u^4-3*u^5);
    yd  = ( y0*(-30*u^2+60*u^3-30*u^4) + y1*(30*u^2-60*u^3+30*u^4) ...
        + T*v0*(1-18*u^2+32*u^3-15*u^4) + T*v1*(-12*u^2+28*u^3-15*u^4) )/T;
    ydd = ( y0*(-60*u+180*u^2-120*u^3) + y1*(60*u-180*u^2+120*u^3) ...
        + T*v0*(-36*u+96*u^2-60*u^3) + T*v1*(-24*u+84*u^2-60*u^3) )/T^2;
end

function [q, qd, qdd] = arm_desired(t, m, dl_deg)
% Prescribed shoulder and elbow motion (rad, rad/s, rad/s^2).
% Shoulder: th_start -> th_rel (clockwise), speeds -w0 -> -w1.
% Elbow   : flexed by dl_deg -> 0 (straight at delivery), starts at rest and
%           is still opening at delivery with rate k*dl/Tb.
    [a, ad, add] = hermite5(t, m.Tb, deg2rad(m.th_start), deg2rad(m.th_rel), -m.w0, -m.w1);
    d = deg2rad(dl_deg);
    [b, bd, bdd] = hermite5(t, m.Tb, d, 0, 0, -m.k*d/m.Tb);
    q   = [a,   b];
    qd  = [ad,  bd];
    qdd = [add, bdd];
end

function [M, C, G] = twoLinkDynamicsMatrices_arm(q, qd, p)
% Class-notes two-link arm (Euler-Lagrange). q = [th1; th2], th1 from the
% horizontal, th2 relative.
%   M(1,1) = m1 lc1^2 + I1 + m2 (L1^2 + lc2^2 + 2 L1 lc2 c2) + I2
%   M(1,2) = m2 (lc2^2 + L1 lc2 c2) + I2
%   M(2,2) = m2 lc2^2 + I2
    m1=p.m1; m2=p.m2; L1=p.L1; lc1=p.lc1; lc2=p.lc2; I1=p.I1; I2=p.I2; g=p.g;
    c2 = cos(q(2));  s2 = sin(q(2));

    M11 = m1*lc1^2 + I1 + m2*(L1^2 + lc2^2 + 2*L1*lc2*c2) + I2;
    M12 = m2*(lc2^2 + L1*lc2*c2) + I2;
    M22 = m2*lc2^2 + I2;
    M = [M11 M12; M12 M22];

    h = -m2*L1*lc2*s2;
    C = [ h*qd(2),   h*(qd(1)+qd(2));
         -h*qd(1),   0              ];

    G1 = (m1*lc1 + m2*L1)*g*cos(q(1)) + m2*lc2*g*cos(q(1)+q(2));
    G2 = m2*lc2*g*cos(q(1)+q(2));
    G = [G1; G2];
end

function tau = arm_inverse_dynamics(t, m, dl_deg, p)
    [q, qd, qdd] = arm_desired(t, m, dl_deg);
    [M,C,G] = twoLinkDynamicsMatrices_arm(q(:), qd(:), p);
    tau = M*qdd(:) + C*qd(:) + G;
end

function dx = twoLinkEOM_arm(t, x, tauFun, p)
    q = x(1:2);  qd = x(3:4);
    [M,C,G] = twoLinkDynamicsMatrices_arm(q, qd, p);
    tau = tauFun(t);
    dx = [qd; M \ (tau - C*qd - G)];
end

function [elbow, tip, vtip] = armKinematics(q, qd, p)
% Positions of elbow and fingertip (shoulder at origin) and fingertip velocity.
    c1 = cos(q(1));  s1 = sin(q(1));
    c12 = cos(q(1)+q(2));  s12 = sin(q(1)+q(2));
    w12 = qd(1) + qd(2);
    elbow = [p.L1*c1; p.L1*s1];
    tip   = elbow + [p.L2*c12; p.L2*s12];
    vtip  = [-p.L1*s1*qd(1) - p.L2*s12*w12;
              p.L1*c1*qd(1) + p.L2*c12*w12];
end

function animateArms(t, R, col, lbl, p, slow)
    nC = numel(R);  N = numel(t);  dtt = t(2) - t(1);
    reach = p.L1 + p.L2;
    lw = linspace(9, 2.5, nC);         % thick (straight arm) -> thin (20 deg), so overlaps stay visible

    fig = figure('Name','Bowling arm animation','Color','w');
    ax = axes(fig); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    xlim(ax, [-1.2 1.2]*reach); ylim(ax, [-1.2 1.2]*reach);
    xlabel(ax,'x (m, bowling direction)'); ylabel(ax,'y (m)');
    plot(ax, 0, 0, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 8);

    h = gobjects(nC,1);  trail = gobjects(nC,1);  ball = gobjects(nC,1);
    for c = 1:nC
        trail(c) = plot(ax, NaN, NaN, '-', 'Color', col(c,:), 'LineWidth', 1);
        h(c)     = plot(ax, [0 0 0], [0 0 0], '-', 'Color', col(c,:), 'LineWidth', lw(c));
        ball(c)  = plot(ax, NaN, NaN, 'o', 'MarkerSize', 8, ...
                        'MarkerFaceColor', col(c,:), 'MarkerEdgeColor', 'k');
    end
    legend(ax, h, lbl, 'Location','southoutside', 'NumColumns', 2);
    titleH = title(ax, '');
    txt = text(ax, -1.15*reach, 1.15*reach, '', 'VerticalAlignment','top', ...
               'FontName','FixedWidth', 'FontSize', 9);

    figure(fig); drawnow; pause(1);          % let the window appear before starting

    for rep = 1:3                            % play the motion three times
        clk = tic;
        while true
            if ~isvalid(fig), return; end
            el = toc(clk)/slow;
            done = el >= t(end);
            if done, k = N; else, k = floor(el/dtt) + 1; end

            s = cell(nC,1);
            for c = 1:nC
                set(h(c), 'XData', [0 R(c).elbow(1,k) R(c).tip(1,k)], ...
                          'YData', [0 R(c).elbow(2,k) R(c).tip(2,k)]);
                set(trail(c), 'XData', R(c).tip(1,1:k), 'YData', R(c).tip(2,1:k));
                set(ball(c),  'XData', R(c).tip(1,k),   'YData', R(c).tip(2,k));
                s{c} = sprintf('%-20s %5.1f m/s', lbl{c}, R(c).speed(k));
            end
            set(txt, 'String', s);
            set(titleH, 'String', sprintf('t = %.3f s   (slow motion x%d)   run %d/3', t(k), slow, rep));
            drawnow limitrate;
            if done, break; end
        end
        pause(0.7);
    end
end