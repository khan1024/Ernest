function retval = inoctave ()
    persistent inout;
    if (isempty (inout))
      inout = exist ('OCTAVE_VERSION', 'builtin') ~= 0;
    end
    retval = inout; 
end