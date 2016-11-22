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

add_executable(testQuadIntersect_CXX testQuadIntersect.cxx)
target_link_libraries(testQuadIntersect_CXX sigrid)

add_executable(testLinearSolve_CXX testLinearSolve.cxx)
target_link_libraries(testLinearSolve_CXX sigrid)

add_executable(testTriangulate_CXX testTriangulate.cxx)
target_link_libraries(testTriangulate_CXX sigrid)

add_executable(testQuadOverlap_CXX testQuadOverlap.cxx)
target_link_libraries(testQuadOverlap_CXX sigrid)

add_executable(testConserveInterp2D_CXX testConserveInterp2D.cxx)
target_link_libraries(testConserveInterp2D_CXX sigrid)

add_executable(testConserveInterp2DBig_CXX testConserveInterp2DBig.cxx)
target_link_libraries(testConserveInterp2DBig_CXX sigrid)

add_executable(testFindOverlappingCells2D_CXX testFindOverlappingCells2D.cxx)
target_link_libraries(testFindOverlappingCells2D_CXX sigrid)

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
add_test(testLinearSolve_CXX ${CMAKE_BINARY_DIR}/tests/testLinearSolve_CXX)
add_test(testQuadIntersect_CXX ${CMAKE_BINARY_DIR}/tests/testQuadIntersect_CXX)
add_test(testTriangulate_CXX ${CMAKE_BINARY_DIR}/tests/testTriangulate_CXX)
add_test(testQuadOverlap_CXX ${CMAKE_BINARY_DIR}/tests/testQuadOverlap_CXX)
add_test(testConserveInterp2D_CXX ${CMAKE_BINARY_DIR}/tests/testConserveInterp2D_CXX)
add_test(testConserveInterp2DBig_CXX ${CMAKE_BINARY_DIR}/tests/testConserveInterp2DBig_CXX)
add_test(testFindOverlappingCells2D_CXX ${CMAKE_BINARY_DIR}/tests/testFindOverlappingCells2D_CXX)

