cmake_minimum_required(VERSION 3.14)
project(kbsa VERSION 0.2.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Paths ---
set(KMC_DIR ${CMAKE_SOURCE_DIR}/extern/kmc)
set(KMC_API_DIR ${KMC_DIR}/kmc_api)
set(KMC_CORE_DIR ${KMC_DIR}/kmc_core)
set(ZLIB_DIR ${KMC_DIR}/3rd_party/cloudflare)

# --- zlib (cloudflare optimized) ---
set(ZLIB_SOURCES
    ${ZLIB_DIR}/adler32.c
    ${ZLIB_DIR}/compress.c
    ${ZLIB_DIR}/crc32.c
    ${ZLIB_DIR}/deflate.c
    ${ZLIB_DIR}/gzclose.c
    ${ZLIB_DIR}/gzlib.c
    ${ZLIB_DIR}/gzread.c
    ${ZLIB_DIR}/gzwrite.c
    ${ZLIB_DIR}/infback.c
    ${ZLIB_DIR}/inffast.c
    ${ZLIB_DIR}/inflate.c
    ${ZLIB_DIR}/inftrees.c
    ${ZLIB_DIR}/trees.c
    ${ZLIB_DIR}/uncompr.c
    ${ZLIB_DIR}/zutil.c
)

add_library(kmc_zlib STATIC ${ZLIB_SOURCES})
target_include_directories(kmc_zlib PUBLIC ${ZLIB_DIR})
set_source_files_properties(${ZLIB_DIR}/deflate.c PROPERTIES COMPILE_FLAGS "-msse4.2")
set_source_files_properties(${ZLIB_DIR}/adler32_simd.c PROPERTIES COMPILE_FLAGS "-mssse3")
set_source_files_properties(${ZLIB_DIR}/crc32_simd.c PROPERTIES COMPILE_FLAGS "-msse4.2 -mpclmul")

# --- KMC API (database reading) ---
set(KMC_API_SOURCES
    ${KMC_API_DIR}/kmc_file.cpp
    ${KMC_API_DIR}/kmer_api.cpp
    ${KMC_API_DIR}/mmer.cpp
)

# --- KMC Core (counting engine) ---
set(KMC_CORE_SOURCES
    ${KMC_CORE_DIR}/bkb_reader.cpp
    ${KMC_CORE_DIR}/bkb_writer.cpp
    ${KMC_CORE_DIR}/cpu_info.cpp
    ${KMC_CORE_DIR}/develop.cpp
    ${KMC_CORE_DIR}/fastq_reader.cpp
    ${KMC_CORE_DIR}/kb_collector.cpp
    ${KMC_CORE_DIR}/kb_completer.cpp
    ${KMC_CORE_DIR}/kb_storer.cpp
    ${KMC_CORE_DIR}/kff_writer.cpp
    ${KMC_CORE_DIR}/kmc_runner.cpp
    ${KMC_CORE_DIR}/kmer.cpp
    ${KMC_CORE_DIR}/mem_disk_file.cpp
    ${KMC_CORE_DIR}/splitter.cpp
    ${KMC_CORE_DIR}/timer.cpp
    ${KMC_CORE_DIR}/rev_byte.cpp
)

# Detect CPU architecture for radix sort
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" HAS_AVX2)
check_cxx_compiler_flag("-mavx" HAS_AVX)
check_cxx_compiler_flag("-msse4.1" HAS_SSE41)

if(HAS_AVX2)
    set_source_files_properties(${KMC_CORE_DIR}/raduls_avx2.cpp PROPERTIES COMPILE_FLAGS "-mavx2")
    list(APPEND KMC_CORE_SOURCES ${KMC_CORE_DIR}/raduls_avx2.cpp)
endif()
if(HAS_AVX)
    set_source_files_properties(${KMC_CORE_DIR}/raduls_avx.cpp PROPERTIES COMPILE_FLAGS "-mavx")
    list(APPEND KMC_CORE_SOURCES ${KMC_CORE_DIR}/raduls_avx.cpp)
endif()
if(HAS_SSE41)
    set_source_files_properties(${KMC_CORE_DIR}/raduls_sse41.cpp PROPERTIES COMPILE_FLAGS "-msse4.1")
    list(APPEND KMC_CORE_SOURCES ${KMC_CORE_DIR}/raduls_sse41.cpp)
endif()
# SSE2 is baseline on x86_64
list(APPEND KMC_CORE_SOURCES ${KMC_CORE_DIR}/raduls_sse2.cpp)

add_library(kmc_core STATIC ${KMC_CORE_SOURCES})
target_include_directories(kmc_core PUBLIC ${KMC_CORE_DIR} ${KMC_DIR} ${ZLIB_DIR})
target_compile_definitions(kmc_core PRIVATE DEVELOP_MODE)
target_link_libraries(kmc_core PUBLIC kmc_zlib)

# --- kbsa executable ---
add_executable(kbsa
    src/main.cpp
    ${KMC_API_SOURCES}
)

target_include_directories(kbsa PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${KMC_API_DIR}
    ${KMC_CORE_DIR}
    ${ZLIB_DIR}
)

target_link_libraries(kbsa PRIVATE kmc_core)

if(UNIX)
    target_link_libraries(kbsa PRIVATE pthread)
endif()

# OpenMP for chromosome-parallel anchor scan
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(kbsa PRIVATE OpenMP::OpenMP_CXX)
endif()

# Release build by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

install(TARGETS kbsa DESTINATION bin)

# --- kbsa_validate (standalone, no KMC dependency) ---
add_executable(kbsa_validate src/kbsa_validate.cpp)
target_compile_features(kbsa_validate PRIVATE cxx_std_17)
install(TARGETS kbsa_validate DESTINATION bin)
