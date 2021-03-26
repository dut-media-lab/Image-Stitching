function compile_mex(debug_flag)
if strfind(computer(),'64')
    defs = '-DA64BITS '; % for 64bit machines - define pointer type 
else      defs = '';
end
if verLessThan('matlab','7.3')          
    defs = [defs, '-DmwIndex=int -DmwSize=size_t ']; 
end
if nargin>0 && debug_flag     
    debugs = ' -g ';     
    sprintf('mex -g'); 
else      debugs = ' -O ';
end
cmd = sprintf('mex %s -largeArrayDims %s texture_mapping.cpp ', debugs, defs);
eval(cmd);

