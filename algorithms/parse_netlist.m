% vim: sw=2:ts=2
function sys=parse_netlist()
  ckt=Ckt;
  sys.E=ckt.C;
  sys.A=-ckt.G;
  sys.B=ckt.B;
end
