% dicom2tif
%   SAM 210802, this script imports the DICOM file sent by Mohammad and outports a TIF.

asdff = [];
for index = 1 : 408
asdf_name = ['E:\Linninger\Subj58\Z',num2str(index + 2696)];
asdf_info = dicominfo(asdf_name);
asdf( :, :, asdf_info.InstanceNumber ) = dicomread(asdf_name);
end
mat2tif(asdf,'E:\Linninger\Subj58\Z2967_to_Z3104.tif')