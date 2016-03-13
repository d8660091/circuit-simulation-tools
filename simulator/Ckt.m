classdef Ckt < handle & dynamicprops

  properties
    name;
    content;
    nodes;
    subckt;
    res;
    cap;
    ind;
    vol;
    cur;
    vccs;
    vcvs;
    mind;
    G;
    C;
    b;
    Q;
    poles;
    reduced_model=cell(1);
    helper;
    pointer; % current bound in G C b matrix, for preallocating memory
    ind_index;
  end

  methods

    function obj=Ckt(varargin)
      try
        if isempty(varargin)
          [file_name,path_name]=uigetfile('sp_files/*.sp','Select circuit to start');
          filename_or_struct=strcat(path_name,file_name);
        else
          filename_or_struct=varargin{1}
        end
      catch err
        disp('input file is not correct');
      end
      if isstruct(filename_or_struct)
        if isfield(filename_or_struct,'name')
          obj.name=filename_or_struct.name;
        else
          obj.name='circuit without metadata';
        end
        obj.G=filename_or_struct.G;
        obj.C=filename_or_struct.C;
        obj.b=filename_or_struct.b;
        if isfield(filename_or_struct,'Q')
          obj.Q=filename_or_struct.Q;
        end
        return;
      end
      obj.helper=Helper(obj);
      ckt_text=fileread(filename_or_struct);
      ckt_text=regexprep(ckt_text,'\s*\n','\n'); %remove ending spaces
      ckt_text=regexprep(ckt_text,'\*.*?\n','\n'); %remove commetn

      % get subckt contents
      subckt_contents=regexp(ckt_text,'\.subckt(.*?).ends','tokens');
      for i=1:length(subckt_contents)  %get subckt
        obj.subckt{i}=Subckt(subckt_contents{i});
      end
      obj.replace_subckt_in_subckt();

      obj.name=regexp(ckt_text, '\<[\w =.]*','match','once'); %get first_line as name of circuit
      tmp=regexprep(ckt_text, '\<[\w =.]*','','once'); %remove first_line
      tmp=regexprep(tmp,'.subckt.*?.ends.*?\s',''); %remove subckt definition
      % Replace subckt instance with general elements
      subckt_inst_str=regexp(tmp,'(?<=\n)X[\w. ]*','match');
      for i=1:length(subckt_inst_str)
        regexprep(subckt_inst_str{i},'\s*$',''); % remove ending spaces
        subckt_name=regexp(subckt_inst_str{i},'\w*$','match','once');
        for j=1:length(obj.subckt)
          if strcmp(obj.subckt{j}.name, subckt_name)
            break;
          end
        end
        replace_expr=obj.subckt{j}.convert2plain(subckt_inst_str{i});
        match_expr=subckt_inst_str{i};
        tmp=regexprep(tmp,match_expr,replace_expr);
      end

      obj.content=tmp;

      res_cell=regexp(tmp,'\n\s*([Rr].*?)(?=[\n\r])','match'); %get resistor
      for i=1:length(res_cell)
        res_elements     = regexp(res_cell{i},'\<[\w.-+]*','match');
        obj.res{i}.name  = res_elements{1};
        obj.res{i}.node1 = res_elements{2};
        obj.res{i}.node2 = res_elements{3};
        obj.res{i}.value = res_elements{4};
      end

      cap_cell=regexp(tmp,'\n\s*([Cc].*?)(?=[\n\r])','match'); %get capacitors
      for i=1:length(cap_cell)
        cap_elements     = regexp(cap_cell{i},'\<[\w.-+]*','match');
        obj.cap{i}.name  = cap_elements{1};
        obj.cap{i}.node1 = cap_elements{2};
        obj.cap{i}.node2 = cap_elements{3};
        obj.cap{i}.value = cap_elements{4};
      end

      ind_cell=regexp(tmp,'\n\s*([Ll].*?)(?=[\n\r])','match'); %get inductors
      for i=1:length(ind_cell)
        ind_elements     = regexp(ind_cell{i},'\<[\w.-+]*','match');
        obj.ind{i}.name  = ind_elements{1};
        obj.ind{i}.node1 = ind_elements{2};
        obj.ind{i}.node2 = ind_elements{3};
        obj.ind{i}.value = ind_elements{4};
      end

      vol_cell=regexp(tmp,'\n\s*([Vv].*?)(?=[\n\r])','match'); %get voltage sources
      for i=1:length(vol_cell)
        vol_elements     = regexp(vol_cell{i},'\<[\w.-+]*','match');
        obj.vol{i}.name  = vol_elements{1};
        obj.vol{i}.node1 = vol_elements{2};
        obj.vol{i}.node2 = vol_elements{3};
        obj.vol{i}.value = vol_elements(4:end);
      end

      cur_cell=regexp(tmp,'\n\s*([Ii].*?)(?=[\n\r])','match'); %get current sources
      for i=1:length(cur_cell)
        cur_elements     = regexp(cur_cell{i},'\<[\w.-+]*','match');
        obj.cur{i}.name  = cur_elements{1};
        obj.cur{i}.node1 = cur_elements{2};
        obj.cur{i}.node2 = cur_elements{3};
        obj.cur{i}.value = cur_elements{end};
      end

      vccs_cell=regexp(tmp,'\n\s*([Gg].*?)(?=[\n\r])','match'); %get VCCS
      for i=1:length(vccs_cell)
        vccs_elements     = regexp(vccs_cell{i},'\<[\w.-+]*','match');
        obj.vccs{i}.name  = vccs_elements{1};
        obj.vccs{i}.nd1   = vccs_elements{2};
        obj.vccs{i}.nd2   = vccs_elements{3};
        obj.vccs{i}.ni1   = vccs_elements{4};
        obj.vccs{i}.ni2   = vccs_elements{5};
        obj.vccs{i}.value = vccs_elements{end};
      end

      vcvs_cell=regexp(tmp,'\n\s*([Ee].*?)(?=[\n\r])','match'); %get VCVS
      for i=1:length(vcvs_cell)
        vcvs_elements     = regexp(vcvs_cell{i},'\<[\w.-+]*','match');
        obj.vcvs{i}.name  = vcvs_elements{1};
        obj.vcvs{i}.nd1   = vcvs_elements{2};
        obj.vcvs{i}.nd2   = vcvs_elements{3};
        obj.vcvs{i}.ni1   = vcvs_elements{4};
        obj.vcvs{i}.ni2   = vcvs_elements{5};
        obj.vcvs{i}.value = vcvs_elements{end};
      end

      mind_cell=regexp(tmp,'\n\s*([Kk].*?)(?=[\n\r])','match'); %get VCVS
      for i=1:length(mind_cell)
        mind_elements     = regexp(mind_cell{i},'\<[\w.-+]*','match');
        obj.mind{i}.name  = mind_elements{1};
        obj.mind{i}.elem1 = mind_elements{2};
        obj.mind{i}.elem2 = mind_elements{3};
        obj.mind{i}.value = mind_elements{end};
      end

      % Make nodes
      obj.nodes=cell(1);
      for i=1:length(obj.res)
        obj.nodes{end+1} = obj.res{i}.node1;
        obj.nodes{end+1} = obj.res{i}.node2;
      end

      for i=1:length(obj.ind)
        obj.nodes{end+1} = obj.ind{i}.node1;
        obj.nodes{end+1} = obj.ind{i}.node2;
      end

      for i=1:length(obj.cap)
        obj.nodes{end+1} = obj.cap{i}.node1;
        obj.nodes{end+1} = obj.cap{i}.node2;
      end

      for i=1:length(obj.cur)
        obj.nodes{end+1} = obj.cur{i}.node1;
        obj.nodes{end+1} = obj.cur{i}.node2;
      end

      for i=1:length(obj.vol)
        obj.nodes{end+1} = obj.vol{i}.node1;
        obj.nodes{end+1} = obj.vol{i}.node2;
      end

      for i=1:length(obj.vcvs)
        obj.nodes{end+1} = obj.vcvs{i}.nd1;
        obj.nodes{end+1} = obj.vcvs{i}.nd2;
        obj.nodes{end+1} = obj.vcvs{i}.ni1;
        obj.nodes{end+1} = obj.vcvs{i}.ni2;
      end

      for i=1:length(obj.vccs)
        obj.nodes{end+1} = obj.vccs{i}.nd1;
        obj.nodes{end+1} = obj.vccs{i}.nd2;
        obj.nodes{end+1} = obj.vccs{i}.ni1;
        obj.nodes{end+1} = obj.vccs{i}.ni2;
      end

      obj.nodes{1}='0';
      obj.nodes=unique(obj.nodes);
      obj.GCb();
      obj.add_mind();
    end

    function GCb(obj)
    % generate G C b matrix
      num_nodes=length(obj.nodes)-1; %exclude ground node
      obj.pointer=num_nodes;
      width=num_nodes+length(obj.ind)+length(obj.vol)+length(obj.vcvs);
      obj.G = sparse(width,width);
      obj.C = sparse(width,width);
      obj.b = sparse(width,1);

      for i=1:length(obj.res)
        obj.add_res(i);
      end

      for i=1:length(obj.ind)
        obj.add_ind(i);
      end

      for i=1:length(obj.cap)
        obj.add_cap(i);
      end

      for i=1:length(obj.cur)
        obj.add_cur(i);
      end

      for i=1:length(obj.vol)
        obj.add_vol(i);
      end

      for i=1:length(obj.vcvs)
        obj.add_vcvs(i);
      end

      for i=1:length(obj.vccs)
        obj.add_vccs(i);
      end
    end

    function index=find_node(obj,nodes_name)
      if ~isa(nodes_name,'integer')
        index=find(ismember(obj.nodes,nodes_name))-1;
      else
        index=nodes_name;
      end
    end

    function add_res(obj, i)
      n1=find(ismember(obj.nodes,obj.res{i}.node1))-1;
      n2=find(ismember(obj.nodes,obj.res{i}.node2))-1;
      val=Helper.convert_value(obj.res{i}.value);
      if (n1 ~= 0)
        obj.G(n1,n1) = obj.G(n1,n1) + 1/val;
      end
      if (n2 ~= 0)
        obj.G(n2,n2) = obj.G(n2,n2) + 1/val;
      end
      if (n1 ~= 0) && (n2 ~= 0)
        obj.G(n1,n2) = obj.G(n1,n2) - 1/val;
        obj.G(n2,n1) = obj.G(n2,n1) - 1/val;
      end
    end

    function add_ind(obj, i)
      n1=find(ismember(obj.nodes,obj.ind{i}.node1))-1;
      n2=find(ismember(obj.nodes,obj.ind{i}.node2))-1;
      val=Helper.convert_value(obj.ind{i}.value);
      d = obj.pointer; % current size of the MNA
      xr = d+1; % new row
      obj.nodes{end+1}=sprintf('I_%s',obj.ind{i}.name);
      % Matlab automatically increases the size of a matrix
      % if you use an index that is obj.bigger than the current size.
      obj.G(xr,xr) = 0; % add new row/column
      obj.C(xr,xr) = 0; % add new row/column
      if (n1 ~= 0)
        obj.G(n1,xr) = 1;
        obj.G(xr,n1) = 1;
      end
      if (n2 ~= 0)
        obj.G(n2,xr) = -1;
        obj.G(xr,n2) = -1;
      end
      obj.C(xr,xr)=-val;
      obj.b(xr)=0;
      obj.ind{i}.index=xr;
      obj.pointer = obj.pointer+1; % current size of the MNA
    end

    function add_cap(obj,i)
      n1=find(ismember(obj.nodes,obj.cap{i}.node1))-1;
      n2=find(ismember(obj.nodes,obj.cap{i}.node2))-1;
      val=Helper.convert_value(obj.cap{i}.value);
      if (n1 ~= 0)
        obj.C(n1,n1) = obj.C(n1,n1) + val;
      end
      if (n2 ~= 0)
        obj.C(n2,n2) = obj.C(n2,n2) + val;
      end
      if (n1 ~= 0) && (n2 ~= 0)
        obj.C(n1,n2) = obj.C(n1,n2) - val;
        obj.C(n2,n1) = obj.C(n2,n1) - val;
      end
    end

    function add_cur(obj,i)
      n1=find(ismember(obj.nodes,obj.cur{i}.node1))-1;
      n2=find(ismember(obj.nodes,obj.cur{i}.node2))-1;
      val=Helper.convert_value(obj.cur{i}.value);
      if (n1 ~= 0)
        obj.b(n1) = obj.b(n1) - val;
      end
      if (n2 ~= 0)
        obj.b(n2) = obj.b(n2) + val;
      end
    end

    function add_vol(obj,i)
      n1=find(ismember(obj.nodes,obj.vol{i}.node1))-1;
      n2=find(ismember(obj.nodes,obj.vol{i}.node2))-1;
      xr = obj.pointer+1; % new row
      obj.b(xr) = 0; % add new row
      % +1 obj.because the existance of ground node.
      obj.nodes{end+1}=sprintf('I_vol(%d)',i); % here +1 needed
      % Matlab automatically increases the size of a matrix
      % if you use an index that is obj.bigger than the current size.
      obj.G(xr,xr) = 0; % add new row/column
      obj.C(xr,xr) = 0; % add new row/column
      if (n1 ~= 0)
        obj.G(n1,xr) = 1;
        obj.G(xr,n1) = 1;
      end
      if (n2 ~= 0)
        obj.G(n2,xr) = -1;
        obj.G(xr,n2) = -1;
      end
      ac_value=find(ismember(lower(obj.vol{i}.value),'ac'))+1;
      if isempty(ac_value)
        ac_value=Helper.convert_value(obj.vol{i}.value{1});
      else
        ac_value=Helper.convert_value(obj.vol{i}.value{ac_value});
      end
      obj.b(xr) = ac_value; % voltage source ac default
      obj.pointer=obj.pointer+1;
    end

    function add_vccs(obj,i)
      nd1=find(ismember(obj.nodes,obj.vccs{i}.nd1))-1;
      nd2=find(ismember(obj.nodes,obj.vccs{i}.nd2))-1;
      ni1=find(ismember(obj.nodes,obj.vccs{i}.ni1))-1;
      ni2=find(ismember(obj.nodes,obj.vccs{i}.ni2))-1;
      val=Helper.convert_value(obj.vccs{i}.value);
      if (nd1 ~= 0) && (ni1 ~= 0 )
        obj.G(nd1,ni1) = obj.G(nd1,ni1) + val;
      end
      if (nd1 ~= 0) && (ni2 ~= 0 )
        obj.G(nd1,ni2) = obj.G(nd1,ni2) - val;
      end
      if (nd2 ~= 0) && (ni1 ~= 0 )
        obj.G(nd2,ni1) = obj.G(nd2,ni1) - val;
      end
      if (nd2 ~= 0) && (ni2 ~= 0 )
        obj.G(nd2,ni2) = obj.G(nd2,ni2) + val;
      end
    end

    function add_vcvs(obj,i)
      nd1=find(ismember(obj.nodes,obj.vcvs{i}.nd1))-1;
      nd2=find(ismember(obj.nodes,obj.vcvs{i}.nd2))-1;
      ni1=find(ismember(obj.nodes,obj.vcvs{i}.ni1))-1;
      ni2=find(ismember(obj.nodes,obj.vcvs{i}.ni2))-1;
      val=Helper.convert_value(obj.vcvs{i}.value);
      d = obj.pointer; % current size of the MNA
      xr = d+1; % new row
      obj.b(xr) = 0; % add new row
      obj.nodes{end+1}=sprintf('I_vcvs(%d)',i);
      % Matlab automatically increases the size of a matrix
      % if you use an index that is obj.bigger than the current size.
      obj.G(xr,xr) = 0; % add new row/column
      obj.C(xr,xr) = 0; % add new row/column
      if (nd1 ~= 0)
        obj.G(nd1,xr) = 1;
        obj.G(xr,nd1) = 1;
      end
      if (nd2 ~= 0)
        obj.G(nd2,xr) = -1;
        obj.G(xr,nd2)= -1;
      end
      if (ni1 ~= 0)
        obj.G(xr,ni1)= -val;
      end
      if (ni2 ~= 0)
        obj.G(xr,ni2)= val;
      end
      obj.pointer=obj.pointer+1;
    end

    function add_mind(obj)
      for i=1:length(obj.mind)
        if isempty(obj.ind_index)
          for j=1:length(obj.ind)
            obj.ind_index.(obj.ind{j}.name).value=str2double(obj.ind{j}.value);
          end
        end
        elem1=obj.mind{i}.elem1;
        elem2=obj.mind{i}.elem2;
        K=Helper.convert_value(obj.mind{i}.value);
        L1=obj.ind_index.(elem1).value;
        L2=obj.ind_index.(elem2).value;
        M=K*sqrt(L1*L2);
        n1=sprintf('I_%s',elem1);
        n1=find(ismember(obj.nodes,n1))-1;
        n2=sprintf('I_%s',elem2);
        n2=find(ismember(obj.nodes,n2))-1;
        obj.C(n1,n2)=obj.C(n1,n2)-M;
        obj.C(n2,n1)=obj.C(n2,n1)-M;
      end
    end

    function bode(obj, n, f ) % Frequency response for output index n, at frequencies f
      addr=sprintf('tmp/%s_original_ac.dat', obj.name);
      if ~exist(addr,'file') || strcmp('Yes', questdlg('Do you want overwrite exisitng file?'))
        s=2*pi*1j*f;
        x=zeros(length(s),1);
        sG=sparse(obj.G);
        sC=sparse(obj.C);
        sb=sparse(obj.b);
        for j=1:length(s)
          X=(sG+s(j)*sC)\sb;
          x(j)=X(n);
        end
        original_ckt_ac=[f',x];
        csvwrite(addr,original_ckt_ac);
      end
    end

    function replace_subckt_in_subckt(obj)
      % remove subckt in subckt
      for i=1:length(obj.subckt)
        subckt_inst_str=regexp(obj.subckt{i}.content,'(?<=\n)X[\w. ]*','match');
        if ~isempty(subckt_inst_str)
          for k=1:length(subckt_inst_str)
            regexprep(subckt_inst_str{i},'\s*$',''); % remove ending spaces in the subckt string
            subckt_name=regexp(subckt_inst_str{i},'\w*$','match','once'); % find subckt class by name
            for j=1:length(obj.subckt)
              if strcmp(obj.subckt{j}.name, subckt_name)
                break;
              end
            end
            replace_expr=obj.subckt{j}.convert2plain(subckt_inst_str{k});
            match_expr=subckt_inst_str{k};
            obj.subckt{i}.content=regexprep(obj.subckt{i}.content,match_expr,replace_expr);
          end
        end
      end
    end

  end
end
