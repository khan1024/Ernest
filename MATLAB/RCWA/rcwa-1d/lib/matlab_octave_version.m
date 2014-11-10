in_octave=inoctave;

old_field_version=2;
if in_octave==0
    % you are working in MATLAB
    version_matlab=[version('-release')];
elseif in_octave==1
    version_octave=[OCTAVE_VERSION];
    number_version_octave=str2num([version_octave(1) version_octave(3) version_octave(5)]);
    if number_version_octave<=324
        old_field_version=1;
    else
        old_field_version=2;
    end
end
