awk 'BEGIN{sum = 0}; $13 > 0 {sum++}; END{print sum}' topo_tec_0151.dat

