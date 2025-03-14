% ODE system of the model of pure and lipidic aSyn aggregation

function xdot = asyn_agg_model(t,x,v)

% ------Variables-------%
M = x(1,1);         % M  = the amount of free monomers in solution
A = x(2,1);         % A  = type-A oligomer number concentration
B = x(3,1);         % B  = type-B oligomer number concentration
F_star = x(4,1);    % Fl  = freshly generated fibril number concentration
Mf_star = x(5,1);   % Mfl = the amount of monomers in freshly generated fibril
Fp = x(6,1);        % Fp = pure fibrils fibril number concentration
Mfp = x(7,1);       % Mfp = the amount of monomers in pure fibrils
Fl = x(8,1);        % Fl  = lipidic fibril number concentration
Mfl = x(9,1);       % Mfl = the amount of monomers in lipidic fibrils
Vfb = x(10,1);      % Vfb = the amount of monomer-coated vesicles bounded to fibrils  
Vf = x(11,1);       % Vf = the amount of free vesicles in solution
Vb = x(12,1);       % Vb  = the amount of monomer-coated lipid vesicles in solution

% ------Parameters-------%
n            = v(1);
kn           = v(2); 
kmax_I       = v(3);
Kp           = v(4);
kmax_II_p    = v(5);
Ks           = v(6);
kc1          = v(7); 
kdis         = v(8);
kc2          = v(9);
kl1          = v(10);           
kl2          = v(11);            
y            = v(12);
kon          = v(13);
KD           = v(14);
k_bind_off   = v(15);
k_bind_on    = v(16);
gamma        = v(17);
Ms           = v(18);
kmax_II_l    = v(19);
kmax_II_star = v(20);

Mftot = Mfl + Mfp - Mf_star;

% -------Equations---------%
M_dot = - n * kn * M.^n - n * kmax_I * M.^n./(M.^n + Kp^n) - n * kmax_II_p * M.^n./(M.^n + Ks^n) .* (Mfp - Mf_star) - n * kmax_II_l * M.^n./(M.^n + Ks^n) .* (Mfl-Mf_star) - n * kmax_II_star * M.^n./(M.^n + Ks^n) .* Mf_star + n * kdis * A + Ms/gamma*(- k_bind_on * M * Vf + k_bind_off * Vb) - kl1 * M * Fp - y * kon * Vb * kl2 * M .* Fl./(kl2 * M + y * kon * Vb);
A_dot = kn * M.^n + kmax_I * M.^n./(M.^n + Kp^n) + kmax_II_p * M.^n./(M.^n + Ks^n) .* (Mfp - Mf_star) + kmax_II_l * M.^n./(M.^n + Ks^n) .* (Mfl - Mf_star) + kmax_II_star * M.^n./(M.^n + Ks^n) .* Mf_star - kc1 * A - kdis * A;
B_dot = kc1 * A - kc2 * B;
F_star_dot = kc2 * B - kl1 * M * F_star - kon * Vb * F_star;
Mf_start_dot = n * kc2 * B - kl1 * M * Mf_star - kon * Vb * Mf_star ;
Fp_dot = kc2 * B - kon * Vb * F_star;
Mfp_dot = n * kc2 * B + kl1 * M * Fp - kon * Vb * Mf_star; 
Fl_dot = kc2 * B - kl1 * M * F_star;
Mfl_dot = n * kc2 * B - kl1 * M * Mf_star + y * kon * Vb * kl2 * M .* Fl./(kl2 * M + y * kon * Vb);
Vf_b_dot = kon * Vb * kl2 * M .* Fl./(kl2 * M + y * kon * Vb);
Vf_dot = - k_bind_on * M * Vf + k_bind_off * Vb;
Vb_dot = k_bind_on * M * Vf - k_bind_off * Vb - kon * Vb * kl2 * M .* Fl./(kl2 * M + y * kon * Vb);

% -------Vector to return-------
xdot = [M_dot; A_dot; B_dot; F_star_dot; Mf_start_dot; Fp_dot; Mfp_dot; Fl_dot; Mfl_dot; Vf_b_dot; Vf_dot; Vb_dot];

end