function m=moments(os,q)
    h=size(os.C,1); % height of each moments
    w=size(os.B,2); % width of each moments block
    m=zeros(h,w,q); % allocate memory for m
    AE=os.A\os.E;
    AB=-os.A\os.B;

    tmp=1;
    for i=1:q
        m(:,:,i)=os.C*tmp*AB;
        if i~=q
            tmp=tmp*AE;
        end
    end % for
