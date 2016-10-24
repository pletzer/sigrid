include_directories(${CMAKE_SOURCE_DIR}/cpp)

configure_file(testLinearSolve.py.in testLinearSolve.py)
configure_file(testLinearSolveOverDetermined.py.in testLinearSolveOverDetermined.py)
configure_file(testTetrahedron.py.in testTetrahedron.py)
configure_file(testTriangle.py.in testTriangle.py)
configure_file(testLine.py.in testLine.py)
configure_file(testBoxIterator.py.in testBoxIterator.py)
configure_file(testFindPointInCell1D.py.in testFindPointInCell1D.py)

add_executable(testFindPointInCell1D_CXX testFindPointInCell1D.cxx)
target_link_libraries(testFindPointInCell1D_CXX sigrid)

add_executable(testFindPointInCellPolar_CXX testFindPointInCellPolar.cxx)
target_link_libraries(testFindPointInCellPolar_CXX sigrid)

add_executable(testFindPointInCellPolarXY_CXX testFindPointInCellPolarXY.cxx)
target_link_libraries(testFindPointInCellPolarXY_CXX sigrid)

add_test(testLinearSolve
        ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testLinearSolve.py)
add_test(testLinearSolveOverDetermined
        ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testLinearSolveOverDetermined.py)
add_test(testBoxIterator
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testBoxIterator.py)
add_test(testTetrahedron
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testTetrahedron.py)
add_test(testTriangle
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testTriangle.py)
add_test(testLine
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testLine.py)
add_test(testFindPointInCell1D
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testFindPointInCell1D.py)
add_test(testFindPointInCell1D_CXX ${CMAKE_BINARY_DIR}/tests/testFindPointInCell1D_CXX)
add_test(testFindPointInCellPolar_CXX ${CMAKE_BINARY_DIR}/tests/testFindPointInCellPolar_CXX)
add_test(testFindPointInCellPolarXY_CXX ${CMAKE_BINARY_DIR}/tests/testFindPointInCellPolarXY_CXX)
