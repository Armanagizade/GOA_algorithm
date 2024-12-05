lowerb = 1;
uperb = 10;
npop  = 50;
nvar = 3;
maxiter = 100;
pop = zeros(2,nvar,npop);
result = zeros(1,maxiter+1);
cmax = 1;
cmin = 0.0001;
c = linspace(cmax,cmin,maxiter);
x_tmp = zeros(npop,nvar);
pos = 1;
fit = 2; 
%random pop-------------------------------------------------------------------------------------------------------
pop(pos,:,:) = unifrnd(lowerb,uperb,1,nvar,npop);
for i=1:npop
    pop(fit,1,i) = sphere(pop(pos,:,i)); 
end
[pppp ,in] = sort(pop(fit,1,:));
pop = pop(:,:,in);
best_solution = pop(:,:,1);
result(1) = best_solution(fit,1);
disp(best_solution);
%main loop -------------------------------------------------------------------------------------------------------
for i=1:maxiter     
    for k=1:npop
        for j=1:npop
            if k ~= j 
            x_tmp(j,:) = c(i) * ((uperb - lowerb)/2) * ( S_func(distance(pop(pos,:,k) , pop(pos,:,j))) .* ((pop(pos,:,k) - pop(pos,:,j)) ./ distance(pop(pos,:,k) , pop(pos,:,j))));
            end
        end
        xnew = c(i) * sum(x_tmp) + best_solution(pos);
        if (all(xnew > lowerb)) & (all(xnew < uperb)) 
            fnew = sphere(xnew);
            
            if best_solution(fit,1) > fnew 
                best_solution(fit,1) = fnew;
                best_solution(pos,:) = xnew;
            
            end
        end
    end
    result(i+1) = best_solution(fit,1);
    disp(i)
    disp(result(i));
end
disp(best_solution)
plot([1:maxiter+1],result)
xlabel('iteration')
ylabel('IAE')
%function --------------------------------------------------------------------------------------------------------
function[y] = sphere(x)
    y = sum(x.^2);
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