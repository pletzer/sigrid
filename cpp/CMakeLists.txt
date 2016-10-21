include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(LIB_FILES
  SgLinearSolve.cpp
  SgBoxIterator.cpp
  SgTetrahedron.cpp
  SgTriangle.cpp
  SgLine.cpp
  SgFindPointInCell.cpp
)

set(HEADER_FILES
  SgLinearSolve.h
  SgBoxIterator.h
  SgTetrahedron.h
  SgTriangle.h
  SgLine.h
  SgNdims.h
  SgFindPointInCell.h
)

add_library(sigrid SHARED ${LIB_FILES})
target_link_libraries(sigrid ${LAPACK_LIBRARIES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install library
install(TARGETS sigrid DESTINATION lib)
