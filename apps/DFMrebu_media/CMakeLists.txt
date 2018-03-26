#
# Add the executable
#

ADD_EXECUTABLE(DFMrebu_media DFMrebu_media.cpp)
TARGET_LINK_LIBRARIES(DFMrebu_media zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin DFMrebu_media)
