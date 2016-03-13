classdef Helper < handle & dynamicprops
  properties
    ckt_handle;
  end

  methods

    function obj=Helper(ckt)
      obj.ckt_handle=ckt;
    end

    function text_formated=format( text ) % Remove starting space and end space
      tmp=regexprep(text, '^\s*', ''); % Remove the starting space for the first line
      text_formated=regexprep(tmp, '\s*\n\s*', '\n'); % Remove the starting and ending space for other line
    end


    function tmp=mor( obj,q )
      if ~isempty(obj.ckt_handle.reduced_model{1})
        qq=zeros(length(obj.ckt_handle.reduced_model),1);
        for i=1:length(obj.ckt_handle.reduced_model)
          qq(i)=obj.ckt_handle.reduced_model{i}.q;
        end
        if ismember(q,qq)
          disp 'Reduced model already exist.'
          return;
        end
      end
      G=obj.ckt_handle.G;
      C=obj.ckt_handle.C;
      b=obj.ckt_handle.b;
      G=sparse(G);
      [L,U,P,Q]=lu(G);
      z=L\(P*b);
      y=U\z;
      R=Q*y; % R=inv(G)*b
      z=L\(P*C);
      y=U\z;
      A=-Q*y; % A=-inv(G)*C
      % Finding orthogonal basis of Krylov subspace
      Q=zeros(size(G,1),q);
      H=zeros(q);
      Q(:,1) = R/norm(R);
      for j = 1:q-1
        Z = A*Q(:,j);
        for i = 1:j
          H(i,j) = Q(:,i).'* Z;
          Z = Z - H(i,j)*Q(:,i);
        end
        H(j+1,j) = norm(Z);
        if H(j+1,j)==0;
          break;
        end
        Q(:,j+1) = Z/H(j+1,j);
      end
      tmp=Reduced_model();
      tmp.q=q;
      tmp.Q=Q;
      tmp.G=Q'*G*Q;
      tmp.C=Q'*C*Q;
      tmp.b=Q'*b;
      tmp.original_ckt=obj.ckt_handle;
      tmp.get_poles();
      if isempty(obj.ckt_handle.reduced_model{1})
        obj.ckt_handle.reduced_model{1}=tmp;
      else
        obj.ckt_handle.reduced_model{end+1}=tmp;
      end
    end

    function [x] = bode(obj,i, n, f ) % Frequency response for output index n, at frequencies f, for reduced_model i
      Q=obj.ckt_handle.reduced_model{i}.Q;
      G=obj.ckt_handle.reduced_model{i}.G;
      C=obj.ckt_handle.reduced_model{i}.C;
      b=obj.ckt_handle.reduced_model{i}.b;
      s=2*pi*1j*f;
      x=zeros(length(f),1);
      for j=1:length(s)
        X=(G+s(j)*C)\b;
        X=Q*X;
        x(j)=X(n);
      end
      semilogx(f,abs(x));
    end

    function p = poles(obj)
      for i=0:length(obj.ckt_handle.reduced_model)
        if i==0
          G=obj.ckt_handle.G;
          C=obj.ckt_handle.C;
          A=-inv(G)*C;
          p=1./nonzeros(eig(full(A)));
          [~,I]=sort(real(p),'descend');
          p=p(I);
          obj.ckt_handle.poles=p;
          tmp=zeros(length(p),length(obj.ckt_handle.reduced_model));
          tmp(1:length(p),i+1)=p;
          continue;
        end
        if isempty(obj.ckt_handle.reduced_model{i})
          break;
        end
        G=obj.ckt_handle.reduced_model{i}.G;
        C=obj.ckt_handle.reduced_model{i}.C;
        A=-inv(G)*C;
        p=1./nonzeros(eig(full(A)));
        [~,I]=sort(real(p),'descend');
        p=p(I);
        obj.ckt_handle.reduced_model{i}.poles=p;
        tmp(1:length(p),i+1)=p;
        p=tmp;
      end
    end
  end

  methods( Static )

    function  [FN]=convert_value( string )
      %CONVERTVALUE Summary of this function goes here
      %   Detailed explanation goes here
      string=lower(string);
      Ostring=regexprep(string,'m','e-3');
      Ostring=regexprep(Ostring,'u','e-6');
      Ostring=regexprep(Ostring,'n','e-9');
      Ostring=regexprep(Ostring,'p','e-12');
      Ostring=regexprep(Ostring,'f','e-15');
      Ostring=regexprep(Ostring,'k','e3');
      Ostring=regexprep(Ostring,'x','e6');
      Ostring=regexprep(Ostring,'g','e9');
      Ostring=regexprep(Ostring,'t','e12');
      FN=sscanf(Ostring,'%f');
    end

    function hspice_data_plot( filename )
      hspice_data=importdata(filename);
      hspice_data=str2double(hspice_data.textdata(3:end,1:3));
      frequence=hspice_data(:,1);
      response=complex(hspice_data(:,2),hspice_data(:,3));
      string=sprintf('Hspice response');
      semilogx(frequence,abs(response),'DisplayName',string);
    end

    function pdata=plot_all( q )
      file_pat=sprintf('tmp/*_%d.dat',q); % Load plotting files
      files=dir(file_pat);
      pdata=struct;
      try % Plot original circuit response
        pdata(1).data=load('tmp/original_ckt_response.dat')
      catch
        error('Original data did not exist');
      end
      semilogx(pdata(1).data(:,1),abs(pdata(1).data(:,2)));
      hold all;
      for i=1:length(files)
        addr=strcat('tmp/',files(i).name);
        pdata(i+1).name=addr;
        pdata(i+1).data=csvread(addr);
        semilogx(pdata(i+1).data(:,1),abs(pdata(i+1).data(:,2)));
      end
    end

    function pdata=plot(varargin)
      p=inputParser;
      default=0;
      addParameter(p,'format','magnitude');
      parse(p,varargin{:});
      [file_name,path_name]=uigetfile('tmp/*.dat','MultiSelect','on');
      pdata=cell(1);
      if ~iscell(file_name)
        addr=strcat(path_name,file_name);
        pdata=csvread(addr);
        subplot(2,1,1);
        semilogx(pdata(:,1),abs(pdata(:,2)));
        subplot(2,1,2);
        semilogx(pdata(:,1),angle(pdata(:,2)));
      else
        for i=1:length(file_name)
          addr=strcat(path_name,file_name{i});
          pdata{i}.name=file_name{i};
          pdata{i}.data=csvread(addr);
          subplot(2,1,1);
          semilogx(pdata{i}.data(:,1),abs(pdata{i}.data(:,2)));
          hold all;
          subplot(2,1,2);
          semilogx(pdata{i}.data(:,1),angle(pdata{i}.data(:,2)));
          hold all;
        end
      end
    end

  end

end
