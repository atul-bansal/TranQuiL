function [optimal_p, tau_in_ns, val, locs, optimal_alpha] = inverseNDFT_robust_alphaReturn(h, f, tau_max,alpha_start, alpha_end)
%% Input
% Simulated Channel h at five different frequencies f. I have assumed that I only have
% single LOS path. I expect a single peak in the Multipath profile. The actual time of flight is 1.82 nanoseconds. 
% h = [exp(1i*-2.4689);exp(1i*0.7381); exp(1i*2.87);exp(1i*-2.98); exp(1i*2.6640)];
% h = [2.4736*exp(1i*-2.1644);0.5292*exp(1i*-1.8426); 2.5779*exp(-1i*1.0209);5.8990*exp(1i*-2.1489); 3.5829*exp(-1i*2.9110)];
% h = [exp(1i*-2.6492);exp(1i*1.5671);exp(-1i*2.3210);exp(1i*0.5509);exp(-1i*2.9724)];
% h = [exp(-1i*1.366);exp(-1i*1.1373); exp(-1i*0.3368);exp(-1i*2.2745); exp(-1i*1.7027)];
% % f = [430e6; 450e6; 520e6; 900e6; 950e6];
% 
% % h = [exp(1i*3.0972);exp(-1i*1.7223);exp(1i*2.0866);exp(1i*0.3803);exp(1i*2.9630);exp(1i*1.1977);exp(-1i*2.3147);exp(1i*0.6358);exp(-1i*0.7929);exp(-1i*2.8043);exp(1i*2.3876);exp(-1i*0.7858);exp(-1i*3.0152);exp(1i*0.7260);exp(1i*2.6267)];
% f = [350e6; 400e6; 450e6; 500e6; 550e6; 600e6; 650e6; 700e6; 750e6; 800e6; 850e6; 910e6; 915e6; 920e6; 950e6];
% sample time of flight values in nanoseconds
% inversef = round((1./f)*1e9);
% sym_inversef = sym(inversef);
% lcm_inversef = lcm(sym_inversef);
% max_tau_ns = double(lcm_inversef)*1e-9;
% tau_in_ns = -max_tau_ns/2:0.01:max_tau_ns/2;

tau_in_ns = -tau_max:0.2:tau_max;
tau_in_ns = tau_in_ns * 1e-9;


% Other parameters taken from Chronos
FourierMatrix = exp(-1i*2*pi*f.*tau_in_ns);
epsilon = 0.001;
% gamma = 1/norm(FourierMatrix)^2;
gamma = 1/norm(FourierMatrix,'fro')^2;
max_peaks = 5;

alpha = alpha_start;
while alpha < alpha_end
p_prev = ones(length(tau_in_ns),1);
p_prev = p_prev / norm(p_prev);
converged = 0;
p_new = 0;
% alpha

%% Function
while converged == 0 && ~isnan(norm(p_new - p_prev))
  p_new = sparsify(p_prev+gamma*FourierMatrix'*(h-FourierMatrix*p_prev),gamma*alpha);
  if norm(p_new - p_prev) < epsilon
    converged = 1;
    optimal_p = p_new;
    fprintf('converged\n');
  else
   norm(p_new - p_prev)
   p_prev = p_new;
   fprintf('Not converged\n');
  end
  1;
end
[val, locs] = findpeaks(abs(optimal_p),tau_in_ns);
if length(locs) > 10*max_peaks
    prev_alpha = alpha;
    alpha = alpha + (length(locs)-max_peaks);
%      alpha = alpha + 50;
    continue;
elseif length(locs) > 5*max_peaks
    prev_alpha = alpha;
    alpha = alpha + (length(locs)-max_peaks);
%     alpha = alpha + 30;
    continue;
elseif length(locs) > 2*max_peaks
    prev_alpha = alpha;
    alpha = alpha + 0.5*(length(locs) - max_peaks);
%     alpha = alpha + 5;
    continue;
elseif length(locs) > max_peaks
    prev_alpha = alpha;
    alpha = alpha + 2;
    continue;
elseif isempty(locs)
    if alpha == alpha_start
        alpha = alpha - 200;
        prev_alpha = 10;
    else
        alpha = (prev_alpha+alpha)/2;
    end
    continue;
else
    length(locs);
    optimal_alpha = alpha;
    break;
end
end
% plot(tau_in_ns, abs(optimal_p));