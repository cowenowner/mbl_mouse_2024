function [Spec] = SPEC_IRASA_spectrogram_compact(Spec)
% FROM: https://purr.purdue.edu/publications/1987/1
% Compacts the structure to save space - these can get big.
% It's not really that useful yet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen - wrapper for the amri_ functions from link above.
% 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 1;
new_fq = Spec.freq(1):res:Spec.freq(end);
new_osci = interp1(Spec.freq,Spec.osci,new_fq,'spline');
new_frac = interp1(Spec.freq,Spec.frac,new_fq,'spline');
new_mixd = interp1(Spec.freq,Spec.mixd,new_fq,'spline');
new_plaw = interp1(Spec.freq,Spec.Plaw,new_fq,'spline');
% Now let's replace.
Spec.freq = single(new_fq);
Spec.osci = single(new_osci);
Spec.frac = single(new_frac);
Spec.mixd = single(new_mixd);
Spec.Plaw = single(new_plaw);
