% calculate the frequency response
function x=freq_resp(os,f)
    n_in=size(os.B,2); % number of inputs
    n_out=size(os.C,1); % number of outputs
    x=zeros(n_in,n_out,length(f));
    for i = 1:n_in
        u=zeros(size(os.B,2),1);
        u(i,1)=1;
        s=2*pi*1j*f;
        for k=1:length(f)
%            tmp=LUsol((s(k)*os.E-os.A),(os.B*u)); % x
            tmp=(s(k)*os.E-os.A)\(os.B*u); % x
            X=os.C*tmp;
            x(i,:,k)=X;
        end % for k
    end % for i
end % response
