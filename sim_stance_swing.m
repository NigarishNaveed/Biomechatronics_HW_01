%% sim_Stance_Swing.m
%
% Simple gait-cycle simulation: run the STANCE simulation, then the SWING
% simulation, and repeat for nCycles. Each phase is exactly the code from
% sim_Leg_Stance.m / sim_Leg_Swing.m, wrapped into a function
% (run_stance, run_swing) with the plots and animation removed.
%
% Why this works: each phase starts from its OWN prescribed start state, and
% the prescribed motions match at the switches:
%   toe-off:     hip -15 deg, knee 5 deg   (end of stance = start of swing)
%   heel-strike: hip +30 deg, knee 5 deg   (end of swing = start of stance)
% -------------------------------------------------------------------

clear; clc; close all;

nCycles = 3;
segs = {}; 
t_all = []; hip_all = []; knee_all = [];
tauHip_all = []; tauKnee_all = [];
GRFx_all = []; GRFy_all = [];            % ground reaction force (stance only)
HRx_all  = []; HRy_all  = [];            % hip reaction force   (swing only)
t0 = 0; tb = [];

for c = 1:nCycles

    % ---------------- stance (from sim_Leg_Stance.m) ----------------
    st = run_stance();
    segs{end+1} = struct('r', st, 'phase', 1);
    n = numel(st.t);
    t_all      = [t_all;      t0 + st.t];
    hip_all    = [hip_all;    st.hip];
    knee_all   = [knee_all;   st.knee];
    tauHip_all = [tauHip_all; st.tau_hip];
    tauKnee_all= [tauKnee_all;st.tau_knee];
    GRFx_all   = [GRFx_all;   st.Rx];
    GRFy_all   = [GRFy_all;   st.Ry];
    HRx_all    = [HRx_all;    nan(n,1)];
    HRy_all    = [HRy_all;    nan(n,1)];
    t0 = t0 + st.t(end);  tb(end+1) = t0;

    % ----------------- swing (from sim_Leg_Swing.m) ------------------
    sw = run_swing();
    segs{end+1} = struct('r', sw, 'phase', 2);
    n = numel(sw.t);
    t_all      = [t_all;      t0 + sw.t];
    hip_all    = [hip_all;    sw.hip];
    knee_all   = [knee_all;   sw.knee];
    tauHip_all = [tauHip_all; sw.tau_hip];
    tauKnee_all= [tauKnee_all;sw.tau_knee];
    GRFx_all   = [GRFx_all;   zeros(n,1)];   % foot is off the ground
    GRFy_all   = [GRFy_all;   zeros(n,1)];
    HRx_all    = [HRx_all;    sw.Rx];
    HRy_all    = [HRy_all;    sw.Ry];
    t0 = t0 + sw.t(end);  tb(end+1) = t0;
end

%% Checks ---------------------------------------------------------------
fprintf('Cycles simulated: %d   total time: %.2f s\n', nCycles, t_all(end));
fprintf('Hip angle at the first switch (toe-off):      stance end %.2f deg, swing start %.2f deg\n', ...
    rad2deg(st.hip(end)), rad2deg(sw.hip(1)));
fprintf('Hip angle at the second switch (heel-strike): swing end %.2f deg, stance start %.2f deg\n', ...
    rad2deg(sw.hip(end)), rad2deg(st.hip(1)));
fprintf('Knee angle at the switches: stance end %.2f, swing start %.2f, swing end %.2f, stance start %.2f deg\n', ...
    rad2deg(st.knee(end)), rad2deg(sw.knee(1)), rad2deg(sw.knee(end)), rad2deg(st.knee(1)));

%% Plots ----------------------------------------------------------------
figure('Name','Gait cycle: angles and torques','Color','w');
subplot(2,1,1);
plot(t_all, rad2deg(hip_all), 'LineWidth',1.5); hold on;
plot(t_all, rad2deg(knee_all), 'LineWidth',1.5); grid on; markPhases(tb);
xlabel('time (s)'); ylabel('angle (deg)');
legend('hip angle','knee angle','Location','best');
title('Joint angles (dotted = switch between stance and swing)');

subplot(2,1,2);
plot(t_all, tauHip_all, 'LineWidth',1.5); hold on;
plot(t_all, tauKnee_all, 'LineWidth',1.5); grid on; markPhases(tb);
xlabel('time (s)'); ylabel('torque (N\cdot m)');
legend('\tau_{hip}','\tau_{knee}','Location','best');
title('Required joint torques');

figure('Name','Gait cycle: reaction forces','Color','w');
plot(t_all, GRFx_all, t_all, GRFy_all, 'LineWidth',1.5); hold on;
plot(t_all, HRx_all, '--', t_all, HRy_all, '--', 'LineWidth',1.5); grid on; markPhases(tb);
xlabel('time (s)'); ylabel('force (N)');
legend('GRF R_x (stance)','GRF R_y (stance)', ...
       'hip reaction R_x (swing)','hip reaction R_y (swing)','Location','best');
title('Ground reaction force (stance) and hip reaction force (swing)');

pk.L1 = 0.45; pk.L2 = 0.45;          % thigh and shank lengths
animateGait(segs, pk, 1);


% =========================================================================
%                         PHASE FUNCTIONS
% =========================================================================

function r = run_stance()
% Everything from sim_Leg_Stance.m except the plots and animation.
    p.ms  = 4.5;  p.mt  = 7.0;  p.mH  = 44.0;
    p.Ls  = 0.45; p.Lt  = 0.45;
    p.lcs = 0.19; p.lct = 0.20;
    p.Ics = 0.05; p.Ict = 0.10;
    p.g   = 9.81;
    T_inv = [1 0; 1 1];  T_mat = [1 0; -1 1];

    Tstance = 0.8;
    hipA_A = 30; hipA_B = -15;
    knee_base = 5; knee_peak = 18;      % knee_base 5 (was 3) to match swing
    desired = @(t) leg_stance_desired(t, Tstance, hipA_A, hipA_B, knee_base, knee_peak);

    tauPhysFun = @(t, state) leg_inverse_dynamics_stance(t, desired, T_inv, p);
    [ph0, phd0, ~] = desired(0);
    state0 = [ph0(1); ph0(2); phd0(1); phd0(2)];
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    [T, Y] = ode45(@(t,s) twoLinkEOM_stance(t, s, tauPhysFun, T_mat, p), [0 Tstance], state0, opts);

    N = length(T);
    tau_hip = zeros(N,1); tau_knee = zeros(N,1); Rx = zeros(N,1); Ry = zeros(N,1);
    for k = 1:N
        [ph_d, phd_d, phdd_d] = desired(T(k));
        tau_phys = leg_inverse_dynamics_stance(T(k), desired, T_inv, p);
        tau_hip(k) = tau_phys(1); tau_knee(k) = tau_phys(2);
        [Rx(k), Ry(k)] = groundReactionForce(ph_d(1), ph_d(2), phd_d(1), phd_d(2), phdd_d(1), phdd_d(2), p);
    end

    r.t = T;
    r.hip = Y(:,1) - Y(:,2);            % hip angle = phi1 - phi2
    r.knee = Y(:,2);
    r.tau_hip = tau_hip; r.tau_knee = tau_knee;
    r.Rx = Rx; r.Ry = Ry;
end

function r = run_swing()
    p.m1  = 7.0;  p.m2  = 4.5;
    p.L1  = 0.45; p.L2  = 0.45;
    p.lc1 = 0.20; p.lc2 = 0.19;
    p.I1  = 0.10; p.I2  = 0.05;
    p.g   = 9.81;

    Tswing = 0.4;
    hip_A = -15; hip_B = 30;
    knee_base = 5; knee_peak = 60; knee_peak_frac = 0.40;
    desired = @(t) leg_swing_desired(t, Tswing, hip_A, hip_B, knee_base, knee_peak, knee_peak_frac);

    tauFun = @(t, state) leg_inverse_dynamics_swing(t, desired, p);
    [th0, thd0, ~] = desired(0);
    state0 = [th0(1); th0(2); thd0(1); thd0(2)];
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    [T, Y] = ode45(@(t,s) twoLinkEOM_swing(t, s, tauFun, p), [0 Tswing], state0, opts);

    N = length(T);
    tau_h = zeros(N,1); tau_k = zeros(N,1); Rx = zeros(N,1); Ry = zeros(N,1);
    for k = 1:N
        [th_d, thd_d, thdd_d] = desired(T(k));
        [M,C,G] = twoLinkDynamicsMatrices_swing(th_d(1), th_d(2), thd_d(1), thd_d(2), p);
        tau = M*thdd_d(:) + C*thd_d(:) + G;
        tau_h(k) = tau(1); tau_k(k) = tau(2);
        [Rx(k), Ry(k)] = hipReactionForce(th_d(1), th_d(2), thd_d(1), thd_d(2), thdd_d(1), thdd_d(2), p);
    end

    r.t = T;
    r.hip = Y(:,1);
    r.knee = Y(:,2);
    r.tau_hip = tau_h; r.tau_knee = tau_k;
    r.Rx = Rx; r.Ry = Ry;
end

function markPhases(tb)
    yl = ylim;
    plot([tb; tb], repmat(yl(:), 1, numel(tb)), 'k:', 'HandleVisibility','off');
end

% =========================================================================
%        FUNCTIONS COPIED UNCHANGED FROM THE TWO ORIGINAL SCRIPTS
% =========================================================================

function [s,sd,sdd] = minjerk(t,T)
    tau = min(max(t,0),T)/T;
    s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;
    sd  = (30*tau.^2 - 60*tau.^3 + 30*tau.^4)/T;
    sdd = (60*tau - 180*tau.^2 + 120*tau.^3)/T^2;
end

% ------------------------------ stance -----------------------------------

function [s,sd,sdd] = sinhump(t,T)
    w = pi/T;
    s   = sin(w*t).^2;
    sd  = 2*sin(w*t).*cos(w*t)*w;
    sdd = 2*w^2*cos(2*w*t);
end

function [ph, phd, phdd] = leg_stance_desired(t, Tstance, hipA_A_deg, hipA_B_deg, knee_base_deg, knee_peak_deg)
    A = deg2rad(hipA_A_deg); B = deg2rad(hipA_B_deg);
    [~,wA,wB] = gaitSpeeds();
    [hipA,hipAd,hipAdd] = hermite5(t, Tstance, A, B, wA, wB);

    base = deg2rad(knee_base_deg); peak = deg2rad(knee_peak_deg);
    [s2,sd2,sdd2] = sinhump(t, Tstance);
    kn   = base + (peak-base)*s2;
    knd  = (peak-base)*sd2;
    kndd = (peak-base)*sdd2;

    phi1   = hipA   + kn;
    phi1d  = hipAd  + knd;
    phi1dd = hipAdd + kndd;
    phi2   = kn; phi2d = knd; phi2dd = kndd;

    ph   = [phi1,   phi2];
    phd  = [phi1d,  phi2d];
    phdd = [phi1dd, phi2dd];
end

function [M, C, G] = twoLinkDynamicsMatrices_stance(phi1, phi2, phi1d, phi2d, p)
    ms=p.ms; mt=p.mt; mH=p.mH; Ls=p.Ls; Lt=p.Lt; lcs=p.lcs; lct=p.lct;
    Ics=p.Ics; Ict=p.Ict; g=p.g;

    M11 = Ics+Ict + Ls^2*(mH+mt) + 2*Ls*(Lt*mH+lct*mt)*cos(phi2) + lcs^2*ms + Lt^2*mH + lct^2*mt;
    M12 = -Ict - Ls*(Lt*mH+lct*mt)*cos(phi2) - Lt^2*mH - lct^2*mt;
    M22 = Ict + Lt^2*mH + lct^2*mt;
    M = [M11 M12; M12 M22];

    h = Ls*(Lt*mH+lct*mt)*sin(phi2);
    C = [-h*phi2d,            h*(phi2d-phi1d);
          h*phi1d,             0             ];

    G1 = -g*( lcs*ms*sin(phi1) + (Lt*mH+lct*mt)*sin(phi1-phi2) + Ls*mH*sin(phi1) + Ls*mt*sin(phi1) );
    G2 =  g*(Lt*mH+lct*mt)*sin(phi1-phi2);
    G = [G1; G2];
end

function tau_phys = leg_inverse_dynamics_stance(t, desired, T_inv, p)
    [ph_d, phd_d, phdd_d] = desired(t);
    [M,C,G] = twoLinkDynamicsMatrices_stance(ph_d(1), ph_d(2), phd_d(1), phd_d(2), p);
    Q = M*phdd_d(:) + C*phd_d(:) + G;
    tau_phys = T_inv * Q;
end

function dstate = twoLinkEOM_stance(t, state, tauPhysFun, T_mat, p)
    phi1 = state(1); phi2 = state(2); phi1d = state(3); phi2d = state(4);
    [M,C,G] = twoLinkDynamicsMatrices_stance(phi1, phi2, phi1d, phi2d, p);
    tau_phys = tauPhysFun(t, state);
    Q = T_mat * tau_phys;
    phidot  = [phi1d; phi2d];
    phiddot = M \ (Q - C*phidot - G);
    dstate = [phidot; phiddot];
end

function [Rx, Ry] = groundReactionForce(phi1, phi2, phi1d, phi2d, phi1dd, phi2dd, p)
    ms=p.ms; mt=p.mt; mH=p.mH; Ls=p.Ls; Lt=p.Lt; lcs=p.lcs; lct=p.lct; g=p.g;

    dphi   = phi1 - phi2;
    ddphi  = phi1d - phi2d;
    dddphi = phi1dd - phi2dd;

    shankTermX = lcs*ms*(phi1dd*cos(phi1) - phi1d^2*sin(phi1));
    thighTermX = mt*( Ls*(phi1dd*cos(phi1)-phi1d^2*sin(phi1)) + lct*(dddphi*cos(dphi)-ddphi^2*sin(dphi)) );
    hatTermX   = mH*( Ls*(phi1dd*cos(phi1)-phi1d^2*sin(phi1)) + Lt*(dddphi*cos(dphi)-ddphi^2*sin(dphi)) );
    Rx = shankTermX + thighTermX + hatTermX;

    shankTermY = lcs*ms*(phi1dd*sin(phi1) + phi1d^2*cos(phi1));
    thighTermY = mt*( Ls*(phi1dd*sin(phi1)+phi1d^2*cos(phi1)) + lct*(dddphi*sin(dphi)+ddphi^2*cos(dphi)) );
    hatTermY   = mH*( Ls*(phi1dd*sin(phi1)+phi1d^2*cos(phi1)) + Lt*(dddphi*sin(dphi)+ddphi^2*cos(dphi)) );
    Ry = (ms+mt+mH)*g - shankTermY - thighTermY - hatTermY;
end

% ------------------------------- swing -----------------------------------

function [th, thd, thdd] = leg_swing_desired(t, Tswing, hip_A_deg, hip_B_deg, knee_base_deg, knee_peak_deg, peak_frac)
    A = deg2rad(hip_A_deg); B = deg2rad(hip_B_deg);
    [~,wA,wB] = gaitSpeeds();
    [th_h,thd_h,thdd_h] = hermite5(t, Tswing, A, B, wB, wA); 

    base = deg2rad(knee_base_deg); peak = deg2rad(knee_peak_deg);
    Tpeak = peak_frac*Tswing;
    if t <= Tpeak
        [s2,sd2,sdd2] = minjerk(t, Tpeak);
        th_k = base + (peak-base)*s2; thd_k = (peak-base)*sd2; thdd_k = (peak-base)*sdd2;
    else
        [s2,sd2,sdd2] = minjerk(t-Tpeak, Tswing-Tpeak);
        th_k = peak + (base-peak)*s2; thd_k = (base-peak)*sd2; thdd_k = (base-peak)*sdd2;
    end

    th   = [th_h,   th_k];
    thd  = [thd_h,  thd_k];
    thdd = [thdd_h, thdd_k];
end

function [M, C, G] = twoLinkDynamicsMatrices_swing(theta_h, theta_k, thetadot_h, thetadot_k, p)
    m1=p.m1; m2=p.m2; L1=p.L1; lc1=p.lc1; lc2=p.lc2; I1=p.I1; I2=p.I2; g=p.g;

    M11 = I1+I2+L1^2*m2+2*L1*lc2*m2*cos(theta_k)+lc1^2*m1+lc2^2*m2;
    M12 = -I2-L1*lc2*m2*cos(theta_k)-lc2^2*m2;
    M22 = I2+lc2^2*m2;
    M = [M11 M12; M12 M22];

    h = L1*lc2*m2*sin(theta_k);
    C = [-h*thetadot_k,               h*(thetadot_k-thetadot_h);
          h*thetadot_h,                0                       ];

    G1 = g*((L1*m2+lc1*m1)*sin(theta_h) + lc2*m2*sin(theta_h-theta_k));
    G2 = -g*lc2*m2*sin(theta_h-theta_k);
    G = [G1; G2];
end

function tau = leg_inverse_dynamics_swing(t, desired, p)
    [th_d, thd_d, thdd_d] = desired(t);
    [M,C,G] = twoLinkDynamicsMatrices_swing(th_d(1), th_d(2), thd_d(1), thd_d(2), p);
    tau = M*thdd_d(:) + C*thd_d(:) + G;
end

function dstate = twoLinkEOM_swing(t, state, tauFun, p)
    theta_h = state(1); theta_k = state(2);
    thetadot_h = state(3); thetadot_k = state(4);
    [M,C,G] = twoLinkDynamicsMatrices_swing(theta_h, theta_k, thetadot_h, thetadot_k, p);
    tau = tauFun(t, state);
    thetadot  = [thetadot_h; thetadot_k];
    thetaddot = M \ (tau - C*thetadot - G);
    dstate = [thetadot; thetaddot];
end

function [Rx, Ry] = hipReactionForce(theta_h, theta_k, thetadot_h, thetadot_k, thetaddot_h, thetaddot_k, p)
    m1=p.m1; m2=p.m2; L1=p.L1; lc1=p.lc1; lc2=p.lc2; g=p.g;

    dphi   = theta_h - theta_k;
    ddphi  = thetadot_h - thetadot_k;
    dddphi = thetaddot_h - thetaddot_k;

    Rx = lc1*m1*(thetaddot_h*cos(theta_h) - thetadot_h^2*sin(theta_h)) ...
       + m2*( L1*(thetaddot_h*cos(theta_h) - thetadot_h^2*sin(theta_h)) ...
            + lc2*( dddphi*cos(dphi) - ddphi^2*sin(dphi) ) );

    Ry = (m1+m2)*g ...
       + lc1*m1*(thetaddot_h*sin(theta_h) + thetadot_h^2*cos(theta_h)) ...
       + m2*( L1*(thetaddot_h*sin(theta_h) + thetadot_h^2*cos(theta_h)) ...
            + lc2*( dddphi*sin(dphi) + ddphi^2*cos(dphi) ) );
end

function [knee, foot] = legFK(hip, kneeAng, p)
% Same kinematics as legForwardKinematics_swing: hip at origin,
% x anterior (+), y up (+).
d1 = [sin(hip), -cos(hip)];
ck = cos(kneeAng); sk = sin(kneeAng);
d2 = [ck*d1(1)+sk*d1(2), -sk*d1(1)+ck*d1(2)];
knee = p.L1*d1;
foot = knee + p.L2*d2;
end

function animateGait(segs, p, slow)
    fps = 60;  reach = p.L1 + p.L2;
    fmax = 0;
    for i = 1:numel(segs), fmax = max(fmax, max(hypot(segs{i}.r.Rx, segs{i}.r.Ry))); end
    forceScale = 0.3*reach/fmax;

    % speed and heel-strike hip height from the first stance segment
    r1 = segs{1}.r;
    [~,f0] = legFK(r1.hip(1),   r1.knee(1),   p);
    [~,f1] = legFK(r1.hip(end), r1.knee(end), p);
    v = (f0(1)-f1(1))/(r1.t(end)-r1.t(1));
    hipHeel = -f0(2);

    fig = figure('Name','Gait Animation','Color','w');
    ax = axes(fig); hold(ax,'on'); axis(ax,'equal'); ylim(ax,[-0.1 1.3]);
    plot(ax,[-50 200],[0 0],'k-','LineWidth',2);
    tx = 0:0.5:20;  plot(ax,tx,zeros(size(tx)),'k|','MarkerSize',10);
    linkLine = plot(ax,[0 0 0],[0 0 0],'b-','LineWidth',3);
    jointPts = plot(ax,[0 0 0],[0 0 0],'ko','MarkerFaceColor','k','MarkerSize',6);
    arrow    = quiver(ax,0,0,0,0,0,'r','LineWidth',2,'MaxHeadSize',0.8);
    titleH   = title(ax,'');
    xlabel(ax,'x (m, walking direction)'); ylabel(ax,'height (m)');
    names = {'STANCE','SWING'};  flipGRFx = 1;

    xf = 0;  hipEnd = [0 0];
    for i = 1:numel(segs)
        r = segs{i}.r;  ph = segs{i}.phase;
        tt  = linspace(r.t(1), r.t(end), max(2, round((r.t(end)-r.t(1))*fps)));
        hip = interp1(r.t,r.hip,tt);  kn = interp1(r.t,r.knee,tt);
        Rx  = interp1(r.t,r.Rx,tt);   Ry = interp1(r.t,r.Ry,tt);
        u   = (tt-tt(1))/(tt(end)-tt(1));
        s   = 10*u.^3 - 15*u.^4 + 6*u.^5;

        for k = 1:numel(tt)
            if ~isvalid(fig), return; end
            [kneeR, footR] = legFK(hip(k), kn(k), p);
            if ph == 1      % foot planted at xf, hip moves over it
                hipW = [xf - footR(1), -footR(2)];
            else            % hip keeps moving at v; height blends to heel-strike height
                hipW = [hipEnd(1) + v*(tt(k)-tt(1)), hipEnd(2) + (hipHeel-hipEnd(2))*s(k)];
            end
            kneeW = hipW + kneeR;  footW = hipW + footR;
            set(linkLine,'XData',[hipW(1) kneeW(1) footW(1)],'YData',[hipW(2) kneeW(2) footW(2)]);
            set(jointPts,'XData',[hipW(1) kneeW(1) footW(1)],'YData',[hipW(2) kneeW(2) footW(2)]);
            if ph == 1
                set(arrow,'XData',footW(1),'YData',footW(2), ...
                    'UData',flipGRFx*Rx(k)*forceScale,'VData',Ry(k)*forceScale);
            else
                set(arrow,'XData',hipW(1),'YData',hipW(2), ...
                    'UData',Rx(k)*forceScale,'VData',Ry(k)*forceScale);
            end
            xlim(ax, hipW(1) + [-1 1]*1.1);           % camera follows the hip
            set(titleH,'String',sprintf('%s  |R| = %.0f N',names{ph},hypot(Rx(k),Ry(k))));
            drawnow; pause(slow/fps);
        end
        if ph == 1, hipEnd = hipW;
        else,       xf = hipW(1) + footR(1);          % next heel-strike location
        end
    end
end

function [y,yd,ydd] = hermite5(t,T,y0,y1,v0,v1)
u = min(max(t,0),T)/T;
y   = y0*(1-10*u^3+15*u^4-6*u^5) + y1*(10*u^3-15*u^4+6*u^5) ...
    + T*v0*(u-6*u^3+8*u^4-3*u^5) + T*v1*(-4*u^3+7*u^4-3*u^5);
yd  = ( y0*(-30*u^2+60*u^3-30*u^4) + y1*(30*u^2-60*u^3+30*u^4) ...
    + T*v0*(1-18*u^2+32*u^3-15*u^4) + T*v1*(-12*u^2+28*u^3-15*u^4) )/T;
ydd = ( y0*(-60*u+180*u^2-120*u^3) + y1*(60*u-180*u^2+120*u^3) ...
    + T*v0*(-36*u+96*u^2-60*u^3) + T*v1*(-24*u+84*u^2-60*u^3) )/T^2;
end

function [v,wA,wB] = gaitSpeeds()
% Forward speed and the hip angular rates needed at heel-strike (wA)
% and toe-off (wB) so the hip keeps moving at speed v through the switches.
L = 0.45; Tstance = 0.8; k0 = deg2rad(5);
hs = deg2rad(30); he = deg2rad(-15);
legX = @(h) L*sin(h) + L*sin(h-k0);     % foot x relative to hip
dX   = @(h) L*cos(h) + L*cos(h-k0);
v  = (legX(hs)-legX(he))/Tstance;       % about 0.86 m/s
wA = -v/dX(hs);   wB = -v/dX(he);       % rad/s
end