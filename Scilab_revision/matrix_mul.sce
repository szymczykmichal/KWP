clear;
A={1 3 5; 2 4 5; 6 6 7}
B={2 4 5; 2 1 1; 6 12 7}
C = zeros(3,3)

for i=1:3
    for j=1:3
        sum = 0
        for k=1:3
           sum = sum + A(i,k)*B(k,j) 
        end
        C(i,j) = sum
    end
end

disp(C);
disp(A*B);
