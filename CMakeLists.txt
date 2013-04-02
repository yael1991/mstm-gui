# CMake project file for our t-matrix framework
cmake_minimum_required (VERSION 2.8)
project (ms-tmatrix)
enable_language(Fortran)

find_package(PythonInterp REQUIRED)

#find any non-standard modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
find_package(Matplotlib REQUIRED)
find_package(PyQt4)

#copy all of the front-end scripts to the binary directory
file(GLOB PY_SCRIPTS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
add_custom_target( frontend-scripts ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${PY_SCRIPTS})
foreach(PY_SCRIPT ${PY_SCRIPTS})
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PY_SCRIPT}
		    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${PY_SCRIPT} ${CMAKE_CURRENT_BINARY_DIR}/${PY_SCRIPT}
		    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${PY_SCRIPTS}
  )
endforeach(PY_SCRIPT)

#copy default input data
configure_file(msinput.inp ${CMAKE_CURRENT_BINARY_DIR}/msinput.inp @ONLY)
configure_file(cube27.pos ${CMAKE_CURRENT_BINARY_DIR}/cube27.pos @ONLY)
configure_file(cyl3000fvp5.pos ${CMAKE_CURRENT_BINARY_DIR}/cyl3000fvp5.pos @ONLY)
configure_file(cylslab3000.pos ${CMAKE_CURRENT_BINARY_DIR}/cylslab3000.pos @ONLY)
configure_file(etaGold.txt ${CMAKE_CURRENT_BINARY_DIR}/etaGold.txt @ONLY)
configure_file(etaSilver.txt ${CMAKE_CURRENT_BINARY_DIR}/etaSilver.txt @ONLY)
configure_file(mstm_guiwindow.ui ${CMAKE_CURRENT_BINARY_DIR}/mstm_guiwindow.ui @ONLY)


#compile the FORTRAN code
file(GLOB SRC "mpidefs-serial.f90" "mstm-intrinsics.f90" "mstm-modules-v2.2.f90" "mstm-main-v2.2.f90")
add_executable(ms-tmatrix ${SRC})