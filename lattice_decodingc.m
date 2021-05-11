function dec_point = lattice_decodingc(x, lattice, cosets)
%this function decodes (quantizes) the point x according to the lattice
%"lattice". Lattices supported are Z^n, D_n, D_n dual, A_n, E_7 dual, E_8
%cosets: the cosets for the definition of E7* via 2Z^7
% Details referred to the paper " J.Conway and N.Sloane, Fast quantizing and decoding algorithms for lattice quantizers and codes, IEEE T-IT, vol.28, no.2, pp.227-232, Mar,1982."