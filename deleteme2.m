w = 1;
C = -0.01;
K = 1;

a = tf([1],[1 2*w*C w.^2]);
step(a);
margin(a);
loopmargin(a,tf(1,[1 0]));
M = feedback(a,tf(1,[1 0]));

[m p om] = bode(a*tf(1,[1 0]));
bode(a*tf(1,[1 0]))
phase = @(om)(-atan((w^2-om.^2)./(-2*C*w*om)));
phase = @(om)(-atan((w^2-om.^2)./(-2*C*w*om)) - pi);
hold on
semilogx(om,phase(om)*180/pi,'r')

amp = @(om)(1./sqrt( (-(1/K)*C*w*om.^2).^2 + (1/K)^2*(w^2*om - om.^3).^2 ))
semilogx(om,20*log10(amp(om)),'r')


Wgm = (2*C*w*tan(pi) - sqrt((-2*C*w*tan(pi))^2 + 4*w^2))/-2;
% Wpm = ;

phase = @(om)(-atand((w^2-om.^2)./(-2*C*w*om)))
plot(-180:1:180,-atan(-180:1:180))
phase = @(om)(-atan((w^2-om.^2)./(2*C*w*om)))
figure
plot(-180:1:180,-atan(-180:1:180))
semilogx(om,phase(om)*180/pi,'r')
phase = @(om)(2*pi-atan((w^2-om.^2)./(2*C*w*om)))
semilogx(om,phase(om)*180/pi,'r')
bode(a*tf(1,[1 0]))
hold on
sm
semilogx(om,phase(om)*180/pi,'r')
figure
semilogx(om,phase(om)*180/pi,'r')
(2*C*w*tan(pi) + sqrt((-2*C*w*tan(pi))^2 + 4*w^2))/-2
margin(a*tf(1,[1 0]))
hold all
margin(2*a*tf(1,[1 0]))
amp = @(om)(sqrt( (-K*C*w*om.^2).^2 + K^2*(w^2*om - om.^3)^2 )
amp = @(om)(sqrt( (-K*C*w*om.^2).^2 + K^2*(w^2*om - om.^3)^2 ))

semilogx(om,amp(om),'r')
amp = @(om)(sqrt( (-K*C*w*om.^2).^2 + K^2*(w^2*om - om.^3).^2 ))
semilogx(om,amp(om),'r')
amp = @(om)(1./sqrt( (-K*C*w*om.^2).^2 + K^2*(w^2*om - om.^3).^2 ))
semilogx(om,amp(om),'r')
semilogx(om,20*log10(amp(om)),'r')
bode(K*a*tf(1,[1 0]))
hold all
semilogx(om,20*log10(amp(om)),'r')



phase = @(om)( pi/2 - atan((-om.^3 + w^2*om) ./ (1 - 2*C*w*om.^2)) );
pow = @(om)(om./sqrt( (1-2*C*w*om.^2).^2 + (-om.^3+w^2*om).^2 ));