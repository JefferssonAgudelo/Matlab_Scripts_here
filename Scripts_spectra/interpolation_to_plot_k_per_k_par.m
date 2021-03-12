t=linspace(0,1,100);
x=sin(t);
tt=linspace(0,2,200);
y=cos(tt);

DP_in_flow_time = interp1(DP_time, DP, flow_time, 'linear','extrap');   % Expresses ?DP? As A Function Of ?flow_time?

xx=interp1(t, x, tt, 'linear','extrap');

figure
plot(tt, y)
hold on
plot(tt, xx)
hold off