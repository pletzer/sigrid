include_directories(${CMAKE_SOURCE_DIR}/cpp
                    ${CMAKE_SOURCE_DIR}/utils)

configure_file(testLinearSolve.py.in testLinearSolve.py)
configure_file(testLinearSolveOverDetermined.py.in testLinearSolveOverDetermined.py)
configure_file(testTetrahedron.py.in testTetrahedron.py)
configure_file(testTriangle.py.in testTriangle.py)
configure_file(testLine.py.in testLine.py)
configure_file(testBoxIterator.py.in testBoxIterator.py)
configure_file(testFindPointInCell1D.py.in testFindPointInCell1D.py)
configure_file(testCreateRotatedPole.py.in testCreateRotatedPole.py)
configure_file(testCreatePolar.py.in testCreatePolar.py)
configure_file(testCreateRectangular.py.in testCreateRectangular.py)

add_executable(testFindPointInCell1D_CXX testFindPointInCell1D.cxx)
target_link_libraries(testFindPointInCell1D_CXX sigrid)

add_executable(testFindPointInCellPolar_CXX testFindPointInCellPolar.cxx)
target_link_libraries(testFindPointInCellPolar_CXX sigrid sgutils)

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
target_link_libraries(testConserveInterp2D_CXX sigrid sgutils)

add_executable(testNodeInterp2D_CXX testNodeInterp2D.cxx)
target_link_libraries(testNodeInterp2D_CXX sigrid sgutils)

add_executable(testConserveInterp2DBig_CXX testConserveInterp2DBig.cxx)
target_link_libraries(testConserveInterp2DBig_CXX sigrid sgutils)

add_executable(testConserveInterpDisk_CXX testConserveInterpDisk.cxx)
target_link_libraries(testConserveInterpDisk_CXX sigrid sgutils)

add_executable(testFindOverlappingCells2D_CXX testFindOverlappingCells2D.cxx)
target_link_libraries(testFindOverlappingCells2D_CXX sigrid sgutils)

add_executable(testSideIterator_CXX testSideIterator.cxx)
target_link_libraries(testSideIterator_CXX sigrid)

add_executable(testNodeInterpRotatedPole testNodeInterpRotatedPole.cxx)
target_link_libraries(testNodeInterpRotatedPole sigrid sgutils)

add_executable(testCmdLineArgParser testCmdLineArgParser.cxx)
target_link_libraries(testCmdLineArgParser sgutils)

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
add_test(testFindPointInCellPolar_1 
	     ${CMAKE_BINARY_DIR}/tests/testFindPointInCellPolar_CXX 
	     -t 1.e-10 -m 100 -p "0, 0.707107" -i "3,10" -nr 21 -nt 65)
add_test(testFindPointInCellPolarXY_CXX ${CMAKE_BINARY_DIR}/tests/testFindPointInCellPolarXY_CXX)
add_test(testLinearSolve_CXX ${CMAKE_BINARY_DIR}/tests/testLinearSolve_CXX)
add_test(testQuadIntersect_CXX ${CMAKE_BINARY_DIR}/tests/testQuadIntersect_CXX)
add_test(testTriangulate_CXX ${CMAKE_BINARY_DIR}/tests/testTriangulate_CXX)
add_test(testQuadOverlap_CXX ${CMAKE_BINARY_DIR}/tests/testQuadOverlap_CXX)
add_test(testConserveInterp2D_CXX ${CMAKE_BINARY_DIR}/tests/testConserveInterp2D_CXX)
add_test(testNodeInterp2D_CXX ${CMAKE_BINARY_DIR}/tests/testNodeInterp2D_CXX)
add_test(testConserveInterp2DBig_CXX ${CMAKE_BINARY_DIR}/tests/testConserveInterp2DBig_CXX)
add_test(testConserveInterpDisk_CXX ${CMAKE_BINARY_DIR}/tests/testConserveInterpDisk_CXX)
add_test(testFindOverlappingCells2D_CXX ${CMAKE_BINARY_DIR}/tests/testFindOverlappingCells2D_CXX)
add_test(testSideIterator_CXX ${CMAKE_BINARY_DIR}/tests/testSideIterator_CXX)
add_test(testNodeInterpRotatedPole_1
         ${CMAKE_BINARY_DIR}/tests/testNodeInterpRotatedPole
         --src_nj 41 --src_ni 81 --dst_nj 21 --dst_ni 41 --delta_lat 30.0 --delta_lon 20)
add_test(testCreateRotatedPole_1
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testCreateRotatedPole.py
         --nj 11 --ni 21 --delta_lat 0.0 --delta_lon 0.0)
add_test(testCreatePolar_1
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testCreatePolar.py
         --nj 11 --ni 21 --radius 1.0)
add_test(testCreateRectangular_1
         ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/testCreateRectangular.py
         --nj 11 --ni 21)

