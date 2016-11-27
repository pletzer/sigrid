include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(LIB_FILES
  CmdLineArgParser.cpp
)

set(HEADER_FILES
  CmdLineArgParser.h
)

add_library(sgutils SHARED ${LIB_FILES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install library
install(TARGETS sgutils DESTINATION lib)
