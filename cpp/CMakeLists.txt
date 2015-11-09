include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIB_FILES
  SgLinearSolve.cpp
  SgBoxIterator.cpp
  SgTetrahedron.cpp
)

set(HEADER_FILES
  SgLinearSolve.h
  SgBoxIterator.h
  SgTetrahedron.h
)

add_library(sigrid SHARED ${LIB_FILES})
target_link_libraries(sigrid ${LAPACK_LIBRARIES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install library
install(TARGETS sigrid DESTINATION lib)
