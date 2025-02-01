clc;
clear;

s0 = read_complex_binary('cohTest_gpsch_0_binary');
s1 = read_complex_binary('cohTest_gpsch_1_binary');

s0_5sec = s0(2e6:2e6+625e3-1); 
s1_5sec = s1(2e6:2e6+625e3-1);

% Take conjugate to find the angle rather than division to remove random
% NaN errors
complex_5sec = angle(s1_5sec.*conj(s0_5sec));
% Create complex number from angles and then take mean of them - correct
% way of averaging angles
mean_angle_complex_5sec = angle(mean(exp(1i*complex_5sec)));
s_coh = s1_5sec + s0_5sec*exp(1i*mean_angle_complex_5sec);
figure; hold on; plot(abs(s_coh))
plot(abs(s0_5sec))
plot(abs(s1_5sec))
legend('Coherent Combined', 'Signal 0','Signal 1');