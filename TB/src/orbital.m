%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef orbital < handle
    properties
        Ham_index  % index in the Hamiltonian matrix
        Rel_index  % index in the 11-orbital basis as defined in PRB 92, 205108 (2015)
        Found_centre % Whether a centre has been assigned
        Centre % Projection of orbital centre onto the metal atoms plane
        %Centre_pris % Projection of orbital centre onto the metal atoms plane for pristine structure
        Pcentre % Orbital centre for p-orbitals (only orbitals centred on Chalcogens)
        %Pcentre_pris % Orbital centre for p-orbitals (only orbitals centred on Chalcogens) for pristine structure
	loc_theta % Local rotation from pristine structure
        Intra_n_list % list of Ham_index of in-plane nearest neighbors
        Intra_n_centre % list of centres (Cartesian coordinates) of in-plane nearest neighbors 
        Intra_n_hoppings % list of centres (Cartesian coordinates) of in-plane nearest neighbors 
        Intra_nn_list % list of Ham_index of in-plane next nearest neighbors
        Intra_nn_centre % list of centres (Cartesian coordinates) of in-plane next nearest neighbors
        Intra_nn_hoppings % list of centres (Cartesian coordinates) of in-plane nearest neighbors 
        Intra_nnn_list % list of Ham_index of in-plane second next nearest neighbors
        Intra_nnn_centre % list of centres (Cartesian coordinates) of in-plane second next nearest neighbors
        Intra_nnn_hoppings % list of centres (Cartesian coordinates) of in-plane nearest neighbors 
        Inter_n_list % list of Ham_index of nearest neighbor on different layers
        Inter_n_centre % list of centres (Cartesian coordinates)  of nearest neighbor on different layers
        Inter_nn_list % list of Ham_index of next nearest neighbor on different layers
        Inter_nn_centre % list of centres (Cartesian coordinates) of next nearest neighbor on different layers
	Inter_mx_n_list % list of Ham_index of nearest neighbor on different layers between metal and chalcogen atoms
        Inter_mx_n_centre % list of centres (Cartesian) of nearest neighbour between M and X on different layers
	Layer % Layer on which orbital's centre lie
        Spin % Spin of orbital
        N % Principal quantum number of corresponding atomic-orbital
        Z % Atomic number of associated atom
        Inner % Whether is an inner chalcogen
    end
    properties (Dependent)
        l % Angular quantum number of corresponding atomic-orbital
        m % Magnetic quantum number of corresponding atomic-orbital
        M1sym % Eigenvalue of symmetry operator M1 on this orbital
        M2sym % Eigenvalue of symmetry operator M2 on this orbital
    end

    % Methods to assign values to class properties
    methods
        function set.Inner(obj,value)
            obj.Inner = value;
        end
        function set.Z(obj,value)
            obj.Z = value;
        end
        function set.Spin(obj,value)
            obj.Spin = value;
        end
        function set.N(obj,value)
            obj.N = value;
        end
        function set.Layer(obj,value)
            obj.Layer = value;
        end
        function set.Intra_n_list(obj,array)
            if(isvector(array))
                obj.Intra_n_list = array;
            else
                error('List is not an array')
            end           
        end
        function set.Intra_nn_list(obj,array)
            if(isvector(array))
                obj.Intra_nn_list = array;
            else
                error('List is not an array')
            end           
         end
         function set.Intra_nnn_list(obj,array)
            if(isvector(array))
                obj.Intra_nnn_list = array;
            else
                error('List is not an array')
            end           
         end
         function set.Intra_n_centre(obj,array)
            %if(isvector(array))
                obj.Intra_n_centre = array;
            %else
            %    error('List is not an array')
            %end           
         end
         function set.Intra_nn_centre(obj,array)
            %if(isvector(array))
                obj.Intra_nn_centre = array;
            %else
            %    error('List is not an array')
            %end           
         end
         function set.Intra_nnn_centre(obj,array)
            %if(isvector(array))
                obj.Intra_nnn_centre = array;
            %else
            %    error('List is not an array')
            %end           
         end
        function set.Intra_n_hoppings(obj,array)
            if(isvector(array))
                obj.Intra_n_hoppings = array;
            else
                error('Hoppings is not an array')
            end           
        end
        function set.Intra_nn_hoppings(obj,array)
            if(isvector(array))
                obj.Intra_nn_hoppings = array;
            else
                error('Hoppings is not an array')
            end           
        end
        function set.Intra_nnn_hoppings(obj,array)
            if(isvector(array))
                obj.Intra_nnn_hoppings = array;
            else
                error('Hoppings is not an array')
            end           
        end
         function set.Inter_n_list(obj,array)
            if(isvector(array))
                obj.Inter_n_list = array;
            else
                error('List is not an array')
            end           
        end
        function set.Inter_n_centre(obj,array)
            %if(isvector(array))
                obj.Inter_n_centre = array;
            %else
            %    error('List is not an array')
            %end           
        end
        function set.Inter_nn_list(obj,array)
            if(isvector(array))
                obj.Inter_nn_list = array;
            else
                error('List is not an array')
            end           
        end
        function set.Inter_nn_centre(obj,array)
            %if(isvector(array))
                obj.Inter_nn_centre = array;
            %else
            %    error('List is not an array')
            %end           
        end
        function set.Inter_mx_n_list(obj,array)
            if(isvector(array))
                obj.Inter_mx_n_list = array;
            else
                error('List is not an array')
            end           
        end
        function set.Inter_mx_n_centre(obj,array)
            %if(isvector(array))
                obj.Inter_mx_n_centre = array;
            %else
            %    error('List is not an array')
            %end           
        end
        function set.Found_centre(obj,value)
            if (islogical(value))
                obj.Found_centre = value;
            else
                error('Found_centre must be a logical variable')
            end
        end
        function set.Ham_index(obj,value)
            obj.Ham_index = value;
        end
        function set.Rel_index(obj,value)
            obj.Rel_index = value;
        end
        function set.Centre(obj,array)
            if(isvector(array))
                obj.Centre = array;
            else
                error('Centre is not an array')
            end
        end
      % function set.Centre_pris(obj,array)
      %     if(isvector(array))
      %         obj.Centre_pris = array;
      %     else
      %         error('Centre_pris is not an array')
      %     end
      % end
        function set.Pcentre(obj,array)
            if(isvector(array))
                obj.Pcentre = array;
            else
                error('Centre is not an array')
            end
        end
      % function set.Pcentre_pris(obj,array)
      %     if(isvector(array))
      %         obj.Pcentre_pris = array;
      %     else
      %         error('Centre_pris is not an array')
      %     end
      % end
        function set.loc_theta(obj,value)
            obj.loc_theta = value;
        end
        %function rindex = get.Rel_index(obj)
        %    if(mod(obj.Ham_index,11) ~= 0)
        %        rindex = mod(obj.Ham_index,11);
        %    else
        %        rindex = 11;
        %    end
        %    if (rindex >11 || rindex < 1)
        %        error('Invalid relative index')
        %    end
        %end
        function m1 = get.M1sym(obj)
            rindex = obj.Rel_index;
            if (rindex > 0 && rindex < 6)
                m1 = -1;
            else
                m1 = +1;
            end
        end
        function m2 = get.M2sym(obj)
            rindex = obj.Rel_index;
            if (rindex == 1 || rindex == 4 || rindex == 7 || rindex == 10)
                m2 = -1;
            else
                m2 = +1;
            end
        end
        function l = get.l(obj)
            rindex = obj.Rel_index;
            if (rindex == 1 || rindex == 2 || rindex == 6 || rindex == 7 || rindex == 8)
                l = 2;
            else
                l = 1;
            end
        end
        function m = get.m(obj)
            rindex = obj.Rel_index;
            if (rindex == 3 || rindex == 9)
                m = 3;
            elseif( rindex == 4 || rindex == 10)
                m = 1;
            elseif( rindex == 5 || rindex == 11)
                m = 2;
            end
        end
    end
end
