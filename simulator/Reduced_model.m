classdef Reduced_model < handle

  properties
    q;
    Q; G; C; b;
    poles;
    second_reduced_model;
    original_ckt;
  end

  methods

    function [x] = bode(obj, n, f ) % Frequency response for output index n, at frequencies f
      s=2*pi*1j*f;
      x=zeros(length(f),1);
      for i=1:length(s)
        X=(obj.G+s(i)*obj.C)\obj.b;
        X=obj.Q*X;
        x(i)=X(n);
      end
      abs_addr=sprintf('tmp/%s_reduced_bode_%d.dat', obj.original_ckt.name, obj.q);
      csvwrite(abs_addr,[f',x]);
    end

    function get_poles(obj)
      A=-inv(obj.G)*obj.C;
      p=1./nonzeros(eig(full(A)));
      [~,I]=sort(real(p),'descend');
      p=p(I);
      obj.poles=p;
    end

    function LSQ=second_level_reduce(obj, node_index, f)
      abs_addr=sprintf('tmp/q_%d_reduced_second.dat', obj.q);
      if exist(abs_addr) || strcmp('No', questdlg('Do you want overwrite exisitng file?'))
        return;
      end
      A=-inv(obj.G)*obj.C;
      [V,D]=eig(A);
      lambda=diag(D);
      index=find(real(lambda<0)); % Index for poles and its eigenvectors on the left hand plane
      Q_second=V(:,index); % Second reduction projection matrix Q
      [Q_second,~]=qr(Q_second);
      Q_second=Q_second(:,1:length(index));
      tmp.Q = Q_second;
      tmp.G = Q_second'*obj.G*Q_second;
      tmp.C = Q_second'*obj.C*Q_second;
      tmp.b = Q_second'*obj.b;
      A=-inv(tmp.G)*tmp.C;  % find poles
      p=1./nonzeros(eig(A));
      [~,I]=sort(real(p),'descend');
      p=p(I);
      tmp.poles=p;
      obj.second_reduced_model=tmp;
      Q_1 = obj.Q; % calculate response
      Q_2 = obj.second_reduced_model.Q;
      s=2*pi*1j*f;
      x_second=zeros(length(f),1);
      for i=1:length(s)
        X=(obj.second_reduced_model.G+s(i)*obj.second_reduced_model.C)\obj.second_reduced_model.b;
        X=Q_1*Q_2*X;
        x_second(i)=X(node_index);
      end
      % Calculate least square error
      x_reduced=zeros(length(f),1);
      for i=1:length(f)
        X=(obj.G+s(i)*obj.C)\obj.b;
        X=obj.Q*X;
        x_reduced(i)=X(node_index);
      end
      LSQ=sum((abs(x_reduced)-abs(x_second)).^2);
      LSQ=LSQ/length(f);
      fprintf('second_level_reduce LSQ = %e\n',LSQ);
      string=sprintf('multi-level reduced model for reduced model %d', obj.q);
      semilogx(f,abs(x_second),'DisplayName',string);
    end

    function iterative_least_square()
      asdf
    end

    function [LSQ,XX]=least_square(obj,node_index,f)
      %obj.original_ckt.helper.poles();
      p=obj.poles;
      n_p_positive=length(p(real(p)>0));
      p=complex(-abs(real(p)),imag(p));
      %p(real(p)>0)=[];
      %p=complex(real(p),abs(imag(p)));
      s_sample=2.*pi.*f*1j;
      A=zeros(length(f),length(p));
      for i=1:length(s_sample)
        for j=1:length(p)
          if(imag(p(j))==0)
            A(i,j)=1./(s_sample(i)-p(j));
          else
            if j>1 && p(j)==conj(p(j-1))
              A(i,j)=1j*1./(s_sample(i)-p(j-1))-1j*1./(s_sample(i)-p(j));
            else
              A(i,j)=1./(s_sample(i)-p(j))+1./(s_sample(i)-p(j+1));
            end
          end
        end
      end
      AA=[real(A);imag(A)];
      % Experimantal value, H(s)_i=x_reduced(i)
      x_reduced=zeros(length(f),1);
      for i=1:length(f)
        X=(obj.G+2*pi*1j*f(i)*obj.C)\obj.b;
        X=obj.Q*X;
        x_reduced(i)=X(node_index);
      end
      HH=[real(x_reduced);imag(x_reduced)];
      XX=AA\HH;
      x_least_square=zeros(length(f),1);
      for i=1:length(s_sample)
        x_least_square(i)=A(i,:)*XX;
      end
      magnitude=abs(x_least_square);
      % Calculate least square error
      LSQ=sum((abs(x_reduced)-magnitude).^2);
      LSQ=LSQ/length(f);
      abs_addr=sprintf('tmp/%s_reduced_least_square_%d.dat',obj.original_ckt.name, obj.q);
      csvwrite(abs_addr,[f',x_least_square]);
    end

  end

end
