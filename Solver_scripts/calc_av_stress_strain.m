function [s_av,e_av] = calc_av_stress_strain(Z,s,e)

s_av = trapz(Z,s,2);
e_av = trapz(Z,e,2);