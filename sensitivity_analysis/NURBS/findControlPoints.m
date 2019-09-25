function xc = findControlPoints(x,y,TOL)

% difference bw consecutive points
x_diff = (x(2:end) - x(1:end-1));

% absolute value of the gradient
grad_y = (y(2:end) - y(1:end-1))./x_diff;
mod_x = x_diff.*abs(grad_y);
plot(abs(grad_y))

% Add end point locations
xc(1) = x(1);
xc(2) = x(end);

% Find points of discontinuity
count = 0;
for i=1:length(grad_y)-1
    if (sign(grad_y(i))~=sign(grad_y(i+1)))
        count = count+1;
        xc = [xc x(i+1)];
    end
end


% Integrate  the absolute gradient
for i=1:length(mod_twist)
    int_y(i) = sum(mod_x(1:i));
end

% Initialize number of control points
no_of_cp = 4;
max_int_y = int_y(end);
error = inf;
while error > TOL
    index = zeros(no_of_cp-1,1);
    for i = 1:no_of_cp-1
         [val , index(i)] =  min(abs(int_y - i*max_int_y/no_of_cp));
         
    end
    
    
end




