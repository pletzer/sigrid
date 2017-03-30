include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(LIB_FILES
  SgLinearSolve.cpp
  SgBoxIterator.cpp
  SgTetrahedron.cpp
  SgTriangle.cpp
  SgLine.cpp
  SgFindPointInCell.cpp
  SgQuadIntersect.cpp
  SgQuadLineIntersect.cpp
  SgQuadLineFlows.cpp
  SgTriangulate.cpp
  SgConserveInterp2D.cpp
  SgFlowInterp2D.cpp
  SgNodeInterp2D.cpp
  SgFindOverlappingCells2D.cpp
  SgOctreePoints.cpp
)

set(HEADER_FILES
  SgLinearSolve.h
  SgBoxIterator.h
  SgTetrahedron.h
  SgTriangle.h
  SgLine.h
  SgNdims.h
  SgFindPointInCell.h
  SgQuadIntersect.h
  SgQuadLineIntersect.h
  SgQuadLineFlows.h
  SgTriangulate.h
  SgConserveInterp2D.h
  SgFlowInterp2D.h
  SgNodeInterp2D.h
  SgFindOverlappingCells2D.h
  SgSideIterator.h
  SgOctreePoints.h
)

add_library(sigrid SHARED ${LIB_FILES})
target_link_libraries(sigrid ${LAPACK_LIBRARIES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install library
install(TARGETS sigrid DESTINATION lib)
