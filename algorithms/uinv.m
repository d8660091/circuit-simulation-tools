% select the column vectors closest to the space of V
function [UinV] = uinv(U, V)
    m=size(U,2);
    n=size(V,2);
    close_value=zeros(m,1);
    for i=1:m
        close_value(i)=max(U(:,i)'*V);
    end
    [~,I]=sort(close_value,'descend');
    I=I(1:n);
    UinV=U(:,I);
    disp('the smallest unrelated value is');
    disp(close_value(n));
end
