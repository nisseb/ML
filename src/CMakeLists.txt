# Example CMakeLists.txt file that uses the FindSeqAn.cmake module for
# building a SeqAn-based app.

project (mlibi_hw1)
cmake_minimum_required (VERSION 2.8.2)

# Place binaries into "bin" directory.
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

# Go to "src" subdirectory.
add_subdirectory (src)

### Same File? ###

# Configure SeqAn, enabling features for libbz2 and zlib.
set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package (SeqAn REQUIRED)

# Add include directories, defines, and flags for SeqAn (and its dependencies).
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Build the program and link it against the SeqAn dependency libraries.
add_executable (main main.cpp)
target_link_libraries (main ${SEQAN_LIBRARIES})