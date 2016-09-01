
configure_file(testLinearSolve.py.in testLinearSolve.py)
configure_file(testLinearSolveOverDetermined.py.in testLinearSolveOverDetermined.py)
configure_file(testTetrahedron.py.in testTetrahedron.py)
configure_file(testTriangle.py.in testTriangle.py)
configure_file(testLine.py.in testLine.py)
configure_file(testBoxIterator.py.in testBoxIterator.py)

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

