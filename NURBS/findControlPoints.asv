function xc = findControlPoints(x,y,TOL)
x_diff = (x(2:end) - x(1:end-1));
grad_y = abs((y(2:end) - y(1:end-1))./x_diff);
mod_x = x_diff.*grad_y;
for i=1:length(mod_twist)
    int_y(i) = sum(mod_x(1:i));
end
Nr_cp = 4;

