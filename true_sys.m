function dydt = true_sys(y,T,A,L)
    
    % y = [y(1) y(2) y(3)] = [(Velocity) (flight path ang) (altitude)]
    % T -> thrust   A -> angle of attack (rad)
    m = 60000;  % mass is 60,000 kg
    g = 9.8;    % accel due to grav is 9.8 ms^(-2)
    
    dydt = zeros(3,1); 
    dydt(1) = (1/m)*(T*cosd(A) - (2.7 + 3.08*(1.25 + 4.2*A)^2)*y(1)^2 - m*g*sind(y(2)));
    dydt(2) = (1/(m*y(1)))*(T*sind(A) + L*68.6*(1.25 + 4.2*A)*y(1)^2 - m*g*cosd(y(2)));
    dydt(3) = y(1)*sind(y(2));
end