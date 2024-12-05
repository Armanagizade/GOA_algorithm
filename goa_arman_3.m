lowerb = 0;
uperb = 50;
npop  = 100;
nvar = 3;
maxiter =100;
pop = zeros(2,nvar,npop);
result = zeros(1,maxiter+1);
cmax = 1;
cmin = 0.0001;
sigmak=0.5*(uperb-lowerb);
pmutation = 0.4;
c = linspace(cmax,cmin,maxiter);
x_tmp = zeros(npop,nvar);
pos = 1;
fit = 2; 
result_thrust = zeros(1,maxiter+1);




pop(pos,:,:) = unifrnd(lowerb,uperb,1,nvar,npop);

for i=1:npop
    pop(fit,1,i) = sphere_thrust(pop(pos,:,i)); 
end
%[pppp ,in] = sort(pop(fit,1,:));
%pop = pop(:,:,in);
best_solution = pop(:,:,1);
result_thrust(1) = best_solution(fit,1);
disp(best_solution);
%main loop -------------------------------------------------------------------------------------------------------
for i=1:maxiter     
    for k=1:npop
        for j=1:npop
            if k ~= j 
            x_tmp(j,:) = c(i) * ((uperb - lowerb)/2) * ( S_func(distance(pop(pos,:,k) , pop(pos,:,j))) .* ((pop(pos,:,k) - pop(pos,:,j)) ./ distance(pop(pos,:,k) , pop(pos,:,j))));
            end

        
        end
        xnew = c(i) * sum(x_tmp) + best_solution(pos,:);
        if randn < pmutation
            
            xnew = xnew + sigmak * randn(1,3);
            
            
        end
        
        if ((all(xnew > lowerb)) & (all(xnew < uperb))) 
            fnew = sphere_thrust(xnew);
            
            if best_solution(fit,1) >= fnew
                disp("x")
                best_solution(fit,1) = fnew;
                pop(pos,:,j) = xnew;
                pop(fit,1,j) = fnew;
                best_solution(pos,:) = xnew;
            
            end
        end
    end
    result_thrust(i+1) = best_solution(fit,1);
    disp(i)
    disp(result_thrust(i));
end
disp(best_solution(1,:))
plot([1:maxiter+1],result_thrust,'LineWidth',1.5);
xlabel('iteration');
ylabel('IAE');




tt = 0:1:maxiter;


function[J] = sphere_thrust(x)
    s = fotf('s');

    B1 = 0.7;
    B2 = 0.1;
    B3 = 0.1;
    B4 = 0.1;
    % roll
    %plant = (65*s+4560)/(s^3 + 109*s^2+1023*s+2935);
    % pitch
    %plant = (56.95*s+4391)/(s^3+105*s^2+870*s+4430);
    %yaw
    %plant = 105/(s^2+413*s);
    %Thrust
    plant = 1.63 / (s^2 + 5*s); 
    Kp = x(1);
    Ki = x(2);
    Ld = 1;
    Kd = x(3);
    Mu = 1;

    cont = fracpid(Kp, Ki, Ld, Kd, Mu);
    sys = feedback(plant*cont, 1);  % Closed-loop system
    dt = 0.01;
    t = 0:dt:10;

    % Compute step response info
    info = stepinfo(step(feedback(plant*cont,1)));
    SettlingTime = info.SettlingTime/10;
    RiseTime = info.RiseTime / 10;
    Overshoot = info.Overshoot;

    % Display settling time, rise time, and overshoot
    
    % cost function

    e = 1 - step(feedback(plant*cont,1),t);

    J = B1 * (sum(t'.*abs(e)*dt)) +  B2 * RiseTime + B3 * SettlingTime + B4 * Overshoot  ;

end
function o=S_func(r)
f=0.5;
l=1.5;
o = ones(1,3);
for m=1:3
    o(m)=f*exp(-r(m)/l)-exp(-r(m));  % Eq. (2.3) in the paper
end
end
function d = distance(a,b)
d = ones(1,3);
for n=1:3
d(n)=sqrt((a(n))^2+(b(n))^2);
end
end
