function y = heaviSideTrue(x)
%Standard version of the heaviSide function (x < 0 = 0, x >= 0 = 1)
y = zeros(size(x, 1), size(x, 2));

for i = 1:length(x)
    if x(i) < 0
        y(i) = 0;
    else
        y(i) = 1;
    end
end

end