cmake_minimum_required(VERSION 3.20)
# C is enabled because the OTIlib sensitivity backend compiles OTIlib's C
# sources directly. Without it CMake silently ignores .c files in a target,
# which manifests only as unresolved symbols at link time.
project(warp_hypad VERSION 0.1.0 LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# warp_hypad -- differentiable warp-spacetime analysis
#
# Milestone 1 only: grid, contiguous storage, fixed-size tensors, 4x4 matrix
# algebra, JSON configuration and metadata export. Milestones 2-10 are not
# implemented; see docs/limitations.md.
#
# No third-party dependencies. Everything here is C++17 and the standard
# library, which keeps the custom-scalar seam under our control -- a vendored
# linear algebra library that could not be instantiated with HyPADScalar would
# have to be removed again in Milestone 10.
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

# --- Provenance -------------------------------------------------------------
#
# Captured at configure time and compiled into the binary, so every
# verification record carries the revision that produced it. A source tree with
# no git history yields "unknown", which the gate reports as a warning rather
# than hiding.

find_package(Git QUIET)
set(WARP_HYPAD_GIT_COMMIT "unknown")
set(WARP_HYPAD_GIT_DIRTY 0)

if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE _git_commit
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _git_result)
  if(_git_result EQUAL 0 AND _git_commit)
    set(WARP_HYPAD_GIT_COMMIT "${_git_commit}")
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} status --porcelain
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE _git_status
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(NOT _git_status STREQUAL "")
    set(WARP_HYPAD_GIT_DIRTY 1)
  endif()
endif()

if(CMAKE_CONFIGURATION_TYPES)
  set(_build_type "$<CONFIG>")
else()
  set(_build_type "${CMAKE_BUILD_TYPE}")
endif()

# --- Core library -----------------------------------------------------------

add_library(warp_hypad_core
  src/core/grid.cpp
  src/core/json.cpp
  src/core/config.cpp
  src/core/metadata.cpp
  src/core/verification.cpp
  src/geometry/finite_difference.cpp
  src/geometry/signature.cpp
  src/metrics/metric.cpp
)

target_include_directories(warp_hypad_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_definitions(warp_hypad_core PRIVATE
  WARP_HYPAD_GIT_COMMIT="${WARP_HYPAD_GIT_COMMIT}"
  WARP_HYPAD_GIT_DIRTY=${WARP_HYPAD_GIT_DIRTY}
  WARP_HYPAD_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
  WARP_HYPAD_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
  WARP_HYPAD_BUILD_TYPE="${_build_type}"
  WARP_HYPAD_PLATFORM="${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
)

# Warnings are errors. A milestone gate that passes while the compiler is
# reporting a signed/unsigned mismatch in an indexing expression is not a gate.
if(MSVC)
  target_compile_options(warp_hypad_core PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  # /fp:precise is the MSVC default and is required, not incidental: the
  # bit-exactness check between double and the custom scalar type assumes the
  # compiler does not reassociate floating-point expressions.
  target_compile_options(warp_hypad_core PUBLIC /fp:precise)
else()
  target_compile_options(warp_hypad_core PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion)
  # Same reasoning as /fp:precise above.
  target_compile_options(warp_hypad_core PUBLIC -ffp-contract=off)
endif()

# --- Verification executables ----------------------------------------------

enable_testing()

add_executable(milestone1_verify
  tests/milestone1_verify.cpp
  tests/grid_checks.cpp
  tests/layout_checks.cpp
  tests/algebra_checks.cpp
  tests/contraction_checks.cpp
  tests/genericity_checks.cpp
  tests/config_checks.cpp
)

target_link_libraries(milestone1_verify PRIVATE warp_hypad_core)
target_include_directories(milestone1_verify PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)

if(MSVC)
  target_compile_options(milestone1_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
else()
  target_compile_options(milestone1_verify PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion)
endif()

set(WARP_HYPAD_RESULTS_DIR "${CMAKE_BINARY_DIR}/results/verification/milestone_01")

# Each group is registered separately so a failure names the area that broke.
foreach(_group grid layout algebra contraction genericity config)
  add_test(NAME milestone1.${_group}
           COMMAND milestone1_verify --group ${_group} --out "${WARP_HYPAD_RESULTS_DIR}")
endforeach()

# The gate itself: the full suite, which is also what writes the reports.
add_test(NAME milestone1.gate
         COMMAND milestone1_verify --out "${WARP_HYPAD_RESULTS_DIR}")

# --- Milestone A2: derivative operators -------------------------------------

add_executable(milestone2_verify
  tests/milestone2_verify.cpp
  tests/derivative_checks.cpp
)
target_link_libraries(milestone2_verify PRIVATE warp_hypad_core)
target_include_directories(milestone2_verify PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
if(MSVC)
  target_compile_options(milestone2_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
else()
  target_compile_options(milestone2_verify PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion)
endif()

set(WARP_HYPAD_M2_RESULTS "${CMAKE_BINARY_DIR}/results/verification/milestone_02")
add_test(NAME milestone2.gate
         COMMAND milestone2_verify --out "${WARP_HYPAD_M2_RESULTS}")

# --- Milestone A2T: temporal derivative operators ---------------------------

add_executable(milestone2t_verify
  tests/milestone2t_verify.cpp
  tests/temporal_checks.cpp
)
target_include_directories(milestone2t_verify PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
if(MSVC)
  target_compile_options(milestone2t_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
else()
  target_compile_options(milestone2t_verify PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion)
endif()

set(WARP_HYPAD_M2T_RESULTS "${CMAKE_BINARY_DIR}/results/verification/milestone_02t")
add_test(NAME milestone2t.gate
         COMMAND milestone2t_verify --out "${WARP_HYPAD_M2T_RESULTS}")
set_tests_properties(milestone2t.gate PROPERTIES DEPENDS "milestone2.gate")

# --- Milestone A3: metric generation ----------------------------------------
#
# Links the OTIlib backend when it is enabled, so the parameter-sensitivity
# checks run against the real HyPAD scalar. With the backend off those cases are
# recorded as SKIPPED rather than silently omitted.

add_executable(milestone3_verify
  tests/milestone3_verify.cpp
  tests/metric_checks.cpp
  tests/metric_sensitivity_checks.cpp
)
target_include_directories(milestone3_verify PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
if(MSVC)
  target_compile_options(milestone3_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
else()
  target_compile_options(milestone3_verify PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion)
endif()

set(WARP_HYPAD_M3_RESULTS "${CMAKE_BINARY_DIR}/results/verification/milestone_03")
add_test(NAME milestone3.gate
         COMMAND milestone3_verify --out "${WARP_HYPAD_M3_RESULTS}")
set_tests_properties(milestone3.gate PROPERTIES DEPENDS "milestone2.gate")

# The gate must run last, so the written record reflects a complete run.
set_tests_properties(milestone1.gate PROPERTIES
  DEPENDS "milestone1.grid;milestone1.layout;milestone1.algebra;milestone1.contraction;milestone1.genericity;milestone1.config")

# --- Plotting smoke test ----------------------------------------------------
#
# Every figure must be reproducible from saved output without rerunning the
# solver. Running the plot script in CTest is what keeps that true: it proves
# the script executes against the files the solver actually wrote, and that no
# figure depends on a hand-edited intermediate.

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import matplotlib"
    RESULT_VARIABLE _mpl_missing
    OUTPUT_QUIET ERROR_QUIET)
  if(_mpl_missing EQUAL 0)
    add_test(NAME milestone1.figures
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/python/plot_milestone1.py"
                     "${WARP_HYPAD_RESULTS_DIR}")
    set_tests_properties(milestone1.figures PROPERTIES DEPENDS "milestone1.gate")

    add_test(NAME milestone2.figures
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/python/plot_milestone2.py"
                     "${WARP_HYPAD_M2_RESULTS}")
    set_tests_properties(milestone2.figures PROPERTIES DEPENDS "milestone2.gate")

    add_test(NAME milestone3.figures
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/python/plot_milestone3.py"
                     "${WARP_HYPAD_M3_RESULTS}")
    set_tests_properties(milestone3.figures PROPERTIES DEPENDS "milestone3.gate")

    add_test(NAME milestone2t.figures
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/python/plot_milestone2t.py"
                     "${WARP_HYPAD_M2T_RESULTS}")
    set_tests_properties(milestone2t.figures PROPERTIES DEPENDS "milestone2t.gate")

    message(STATUS "  figures:  enabled (${Python3_EXECUTABLE})")
  else()
    message(STATUS "  figures:  SKIPPED (matplotlib not importable)")
  endif()
else()
  message(STATUS "  figures:  SKIPPED (no Python 3 interpreter found)")
endif()

# --- Program B, Milestone B1: EM effective-metric module --------------------
#
# A separate module with its own namespace and include root. It depends on
# warp_hypad_core for Mat4 and the verification harness, and on nothing else --
# in particular it does NOT depend on Program A's curvature pipeline, which is
# the point: B1 proceeds in parallel and is not blocked by it.

set(PROGRAM_B_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../program_b_em_inverse_design")

if(EXISTS "${PROGRAM_B_ROOT}/effective_metric/src/metric_medium_map.cpp")
  add_library(program_b_em STATIC
    "${PROGRAM_B_ROOT}/effective_metric/src/metric_medium_map.cpp"
    "${PROGRAM_B_ROOT}/effective_metric/src/covariant_constitutive.cpp"
  )
  target_include_directories(program_b_em PUBLIC
    $<BUILD_INTERFACE:${PROGRAM_B_ROOT}/constitutive/include>
    $<BUILD_INTERFACE:${PROGRAM_B_ROOT}/effective_metric/include>
    $<BUILD_INTERFACE:${PROGRAM_B_ROOT}/inverse_design/include>
  )
  target_sources(program_b_em PRIVATE
    "${PROGRAM_B_ROOT}/inverse_design/src/design_strings.cpp"
  )
  target_link_libraries(program_b_em PUBLIC warp_hypad_core)
  if(MSVC)
    target_compile_options(program_b_em PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(program_b_em PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  add_executable(b1_verify
    "${PROGRAM_B_ROOT}/verification/b1_verify.cpp"
    "${PROGRAM_B_ROOT}/verification/b1_checks.cpp"
  )
  target_include_directories(b1_verify PRIVATE "${PROGRAM_B_ROOT}/verification")
  target_link_libraries(b1_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b1_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b1_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  add_executable(b1c_s23_verify
    "${PROGRAM_B_ROOT}/verification/b1c_s23_verify.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_s23_checks.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_s5_checks.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_dual_checks.cpp"
  )
  target_include_directories(b1c_s23_verify PRIVATE "${PROGRAM_B_ROOT}/verification")
  target_link_libraries(b1c_s23_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b1c_s23_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b1c_s23_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # Gate B1C-S456: field reconstruction (step 4), the magnetoelectric
  # remainder (step 5), and Z reconstruction from the observer blocks via
  # (A11) (step 6). Runs after B1C-S23, which it takes as given.
  add_executable(b1c_s456_verify
    "${PROGRAM_B_ROOT}/verification/b1c_s456_verify.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_s456_checks.cpp"
  )
  target_include_directories(b1c_s456_verify PRIVATE "${PROGRAM_B_ROOT}/verification")
  target_link_libraries(b1c_s456_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b1c_s456_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b1c_s456_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  add_executable(b1c_verify
    "${PROGRAM_B_ROOT}/verification/b1c_verify.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_checks.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_velocity_checks.cpp"
    "${PROGRAM_B_ROOT}/verification/b1c_dual_checks.cpp"
  )
  target_include_directories(b1c_verify PRIVATE "${PROGRAM_B_ROOT}/verification")
  target_link_libraries(b1c_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b1c_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b1c_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # --- Program B, Milestone B2: abstract constitutive inverse design ---------
  #
  # Optimizes over the A11 reconstruction verified at B1C-S456. It carries no new
  # physics conventions, so it depends on the B1C gates through CTest rather than
  # repeating their cases.
  add_executable(b2_verify
    "${PROGRAM_B_ROOT}/inverse_design/verification/b2_verify.cpp"
    "${PROGRAM_B_ROOT}/inverse_design/verification/b2_checks.cpp"
  )
  target_include_directories(b2_verify PRIVATE
    "${PROGRAM_B_ROOT}/inverse_design/verification")
  target_link_libraries(b2_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b2_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b2_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # --- Program B, Milestone B3: physical admissibility ----------------------
  add_executable(b3_verify
    "${PROGRAM_B_ROOT}/inverse_design/verification/b3_verify.cpp"
    "${PROGRAM_B_ROOT}/inverse_design/verification/b3_checks.cpp"
  )
  target_include_directories(b3_verify PRIVATE
    "${PROGRAM_B_ROOT}/inverse_design/verification")
  target_link_libraries(b3_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b3_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b3_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # --- Program B, Milestone B4: Maxwell propagation -------------------------
  add_executable(b4_verify
    "${PROGRAM_B_ROOT}/inverse_design/verification/b4_verify.cpp"
    "${PROGRAM_B_ROOT}/inverse_design/verification/b4_checks.cpp"
  )
  target_include_directories(b4_verify PRIVATE
    "${PROGRAM_B_ROOT}/inverse_design/verification")
  target_link_libraries(b4_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b4_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b4_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # --- Program B, Milestone B4C: propagation correction ---------------------
  add_executable(b4c_verify
    "${PROGRAM_B_ROOT}/inverse_design/verification/b4c_verify.cpp"
    "${PROGRAM_B_ROOT}/inverse_design/verification/b4c_checks.cpp"
  )
  target_include_directories(b4c_verify PRIVATE
    "${PROGRAM_B_ROOT}/inverse_design/verification")
  target_link_libraries(b4c_verify PRIVATE program_b_em)
  if(MSVC)
    target_compile_options(b4c_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(b4c_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  set(WARP_HYPAD_B1_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b1")
  set(WARP_HYPAD_B1C_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b1c")
  add_test(NAME programb1.gate COMMAND b1_verify --out "${WARP_HYPAD_B1_RESULTS}")
  set(WARP_HYPAD_B1CS23_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b1c_s23")
  add_test(NAME programb1c_s23.gate
           COMMAND b1c_s23_verify --out "${WARP_HYPAD_B1CS23_RESULTS}")
  set(WARP_HYPAD_B1CS456_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b1c_s456")
  add_test(NAME programb1c_s456.gate
           COMMAND b1c_s456_verify --out "${WARP_HYPAD_B1CS456_RESULTS}")
  set_tests_properties(programb1c_s456.gate PROPERTIES DEPENDS "programb1c_s23.gate")
  add_test(NAME programb1c.gate COMMAND b1c_verify --out "${WARP_HYPAD_B1C_RESULTS}")
  # One call, not three: DEPENDS is a property, so successive set_tests_properties
  # calls REPLACE it rather than accumulate. Two of the three dependencies were
  # being silently dropped.
  set_tests_properties(programb1c.gate PROPERTIES
    DEPENDS "programb1.gate;programb1c_s23.gate;programb1c_s456.gate")

  set(WARP_HYPAD_B2_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b2")
  add_test(NAME programb2.gate COMMAND b2_verify --out "${WARP_HYPAD_B2_RESULTS}")
  # One assignment, listing every dependency: DEPENDS is a property, so
  # successive set_tests_properties calls replace it rather than appending.
  set_tests_properties(programb2.gate PROPERTIES
    DEPENDS "programb1c_s23.gate;programb1c_s456.gate;programb1c.gate")

  set(WARP_HYPAD_B3_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b3")
  add_test(NAME programb3.gate COMMAND b3_verify --out "${WARP_HYPAD_B3_RESULTS}")
  set_tests_properties(programb3.gate PROPERTIES
    DEPENDS "programb1c.gate;programb2.gate")

  set(WARP_HYPAD_B4_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b4")
  add_test(NAME programb4.gate COMMAND b4_verify --out "${WARP_HYPAD_B4_RESULTS}")
  set_tests_properties(programb4.gate PROPERTIES
    DEPENDS "programb2.gate;programb3.gate")

  set(WARP_HYPAD_B4C_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_b4c")
  add_test(NAME programb4c.gate COMMAND b4c_verify --out "${WARP_HYPAD_B4C_RESULTS}")
  set_tests_properties(programb4c.gate PROPERTIES DEPENDS "programb4.gate")
  message(STATUS "  program B: B1 enabled")
else()
  message(STATUS "  program B: not present, skipped")
endif()

# --- V0.1A: Program A swept-field exporter ----------------------------------
#
# The exporter MAY evaluate physics: it calls the same AlcubierreMetric path the
# numerical gates use. The Python renderers may not, and the visual gate greps
# their AST to keep it that way.
add_executable(a_sweep_export tools/a_sweep_export.cpp)
target_link_libraries(a_sweep_export PRIVATE warp_hypad_core)
if(MSVC)
  target_compile_options(a_sweep_export PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
else()
  target_compile_options(a_sweep_export PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
# BUILD-LOCAL outputs. The OTI-enabled and OTI-disabled configurations must not
# write into the same directory: they produce different sensitivity backends, and
# whichever suite ran last would silently define the committed record. Each build
# consumes its own sweeps, and running either suite leaves the source tree clean.
set(WARP_HYPAD_A_SWEEPS
    "${CMAKE_BINARY_DIR}/results/verification/v0_program_a_sweeps")
set(WARP_HYPAD_A_GALLERY
    "${CMAKE_BINARY_DIR}/results/visualization/program_a")
# The canonical committed records. Nothing writes here automatically.
set(WARP_HYPAD_M3_RESULTS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/results/verification/milestone_03")
set(WARP_HYPAD_A_SWEEPS_CANONICAL
    "${CMAKE_CURRENT_SOURCE_DIR}/results/verification/v0_program_a_sweeps")
set(WARP_HYPAD_A_GALLERY_CANONICAL
    "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/galleries/program_a")

# Explicit, opt-in promotion. A test must never do this: promotion is a decision
# about which verified build defines the committed record, and a test that made
# it silently would let the last suite to run rewrite the repository.
#
#   cmake --build build --target promote_a_sweeps
add_custom_target(promote_a_sweeps
  COMMAND ${CMAKE_COMMAND} -E make_directory "${WARP_HYPAD_A_SWEEPS_CANONICAL}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${WARP_HYPAD_A_SWEEPS}" "${WARP_HYPAD_A_SWEEPS_CANONICAL}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${WARP_HYPAD_A_GALLERY_CANONICAL}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          "${WARP_HYPAD_A_GALLERY}" "${WARP_HYPAD_A_GALLERY_CANONICAL}"
  COMMENT "Promoting this build's Program A sweeps and visuals to the source tree"
  VERBATIM)
add_test(NAME programa.visual_sweep_export
         COMMAND a_sweep_export --out "${WARP_HYPAD_A_SWEEPS}")
if(Python3_Interpreter_FOUND)
  add_test(NAME programa.visual_static COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/program_a/plot_metric_fields.py"
    --results "${WARP_HYPAD_M3_RESULTS_SRC}"
    --out "${WARP_HYPAD_A_GALLERY}")
  add_test(NAME programa.visual_render COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/program_a/animate_sweeps.py"
    --sweeps "${WARP_HYPAD_A_SWEEPS}"
    --out "${WARP_HYPAD_A_GALLERY}/animations"
    --contact-sheets "${WARP_HYPAD_A_GALLERY}/contact_sheets")
  add_test(NAME programa.visual_gate COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/regression/v0_program_a_verify.py"
    --gallery "${WARP_HYPAD_A_GALLERY}"
    --sweeps "${WARP_HYPAD_A_SWEEPS}"
    --visualization "${CMAKE_CURRENT_SOURCE_DIR}/../visualization")
  # One DEPENDS assignment per test: repeated set_tests_properties(... DEPENDS)
  # calls REPLACE the property rather than appending to it.
  # One DEPENDS assignment per test: repeated set_tests_properties(... DEPENDS)
  # calls REPLACE the property rather than appending, so the gate names both
  # producers in a single list.
  set_tests_properties(programa.visual_render PROPERTIES
    DEPENDS "programa.visual_sweep_export")
  set_tests_properties(programa.visual_gate PROPERTIES
    DEPENDS "programa.visual_render;programa.visual_static")
endif()

# --- Program C, Milestone C0: catalogue-grounded astrogation -----------------
#
# A SEPARATE program. It gets its own directory, namespace and gate, and Programs
# A and B know nothing about it. The dependency edge runs one way only, and at C0
# program_c consumes only warp_hypad_core for the verification harness -- no
# catalogue, route or UI concern enters an A or B core.
#
# It builds from this file rather than its own because CTest is unified here, the
# same way program_b does. That is a build-tree convenience, not a coupling.

set(PROGRAM_C_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../program_c_astrogation")

if(EXISTS "${PROGRAM_C_ROOT}/verification/c0_verify.cpp")
  add_executable(c0_verify
    "${PROGRAM_C_ROOT}/verification/c0_verify.cpp"
    "${PROGRAM_C_ROOT}/verification/c0_checks.cpp"
  )
  target_include_directories(c0_verify PRIVATE
    "${PROGRAM_C_ROOT}/include"
    "${PROGRAM_C_ROOT}/verification")
  target_link_libraries(c0_verify PRIVATE warp_hypad_core)
  if(MSVC)
    target_compile_options(c0_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(c0_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  add_executable(c0_demo "${PROGRAM_C_ROOT}/demo/c0_demo.cpp")
  target_include_directories(c0_demo PRIVATE "${PROGRAM_C_ROOT}/include"
    "${PROGRAM_C_ROOT}/verification")
  target_link_libraries(c0_demo PRIVATE warp_hypad_core)
  if(MSVC)
    target_compile_options(c0_demo PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(c0_demo PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # The real-data gate. Reads the COMMITTED Gaia DR3 snapshot -- offline and
  # deterministic, like every other gate here. Live retrieval lives only in
  # tools/import_gaia.py, which no test invokes.
  add_executable(c0_real_verify "${PROGRAM_C_ROOT}/verification/c0_real_verify.cpp")
  target_include_directories(c0_real_verify PRIVATE "${PROGRAM_C_ROOT}/include")
  target_link_libraries(c0_real_verify PRIVATE warp_hypad_core)
  if(MSVC)
    target_compile_options(c0_real_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(c0_real_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()

  # --- C1: relativistic astrogation ---------------------------------------
  if(EXISTS "${PROGRAM_C_ROOT}/verification/c1_verify.cpp")
    add_executable(c1_verify "${PROGRAM_C_ROOT}/verification/c1_verify.cpp")
    target_include_directories(c1_verify PRIVATE "${PROGRAM_C_ROOT}/include")
    target_link_libraries(c1_verify PRIVATE warp_hypad_core)
    if(MSVC)
      target_compile_options(c1_verify PRIVATE /W4 /WX /permissive- /Zc:__cplusplus)
    else()
      target_compile_options(c1_verify PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    set(WARP_HYPAD_C1_RESULTS "${PROGRAM_C_ROOT}/results/verification/program_c1")
    add_test(NAME programc1.gate COMMAND c1_verify
      --snapshot "${PROGRAM_C_ROOT}/snapshots/gaia_dr3_10pc.txt"
      --out "${WARP_HYPAD_C1_RESULTS}")
    # The visualization chain, in order: the solver writes the records, the
    # driver draws them, the regression gate checks them. CTest DEPENDS enforces
    # the ordering, so a visual can never be checked against a stale record.
    if(Python3_Interpreter_FOUND)
      add_test(NAME programc1.visuals COMMAND "${Python3_EXECUTABLE}"
        "${PROGRAM_C_ROOT}/visualization/make_c1_visuals.py"
        --results "${WARP_HYPAD_C1_RESULTS}")
      add_test(NAME programc1.visualgate COMMAND "${Python3_EXECUTABLE}"
        "${PROGRAM_C_ROOT}/visualization/tests/c1_visual_verify.py"
        --results "${WARP_HYPAD_C1_RESULTS}"
        --visualization "${PROGRAM_C_ROOT}/visualization")
      # One assignment per test: repeated set_tests_properties(... DEPENDS ...)
      # calls REPLACE the property rather than appending to it.
      set_tests_properties(programc1.visuals PROPERTIES DEPENDS "programc1.gate")
      set_tests_properties(programc1.visualgate PROPERTIES
        DEPENDS "programc1.visuals")
    endif()
    message(STATUS "  program C: C1 enabled")
  endif()

  set(WARP_HYPAD_C0_RESULTS "${CMAKE_BINARY_DIR}/results/verification/program_c0")
  add_test(NAME programc0.gate COMMAND c0_verify --out "${WARP_HYPAD_C0_RESULTS}")
  add_test(NAME programc0.demo COMMAND c0_demo --out "${WARP_HYPAD_C0_RESULTS}")
  set_tests_properties(programc0.demo PROPERTIES DEPENDS "programc0.gate")
  add_test(NAME programc0.real COMMAND c0_real_verify
    --snapshot "${PROGRAM_C_ROOT}/snapshots/gaia_dr3_5pc.txt"
    --out "${CMAKE_BINARY_DIR}/results/verification/program_c0_real")
  message(STATUS "  program C: C0 enabled")
else()
  message(STATUS "  program C: not present, skipped")
endif()

# --- OTIlib sensitivity backend (opt-in) ------------------------------------
#
# OFF by default, and that is an ARCHITECTURAL choice, not an open licensing
# question. GPL-3.0 is accepted for the research prototype; the decision is made.
#
# What opt-in buys is that the core library and every non-sensitivity gate build
# and pass without the dependency, so the suite stays runnable for anyone who
# does not want it, and the OTIlib-disabled configuration is checked on every
# release. The backend is isolated behind include/warp_hypad/sensitivity/ and
# pinned externally rather than vendored.
#
#   cmake -S . -B build -DWARP_HYPAD_WITH_OTILIB=ON
#
# Technical note: upstream states OTIlib on Windows works "only under WSL". That
# applies to the full CMake/conda build with its `make gendata` step. The static
# scalar type onumm3n2 has no precomputed-data dependency, so the scalar-only
# translation unit below builds natively with MSVC. Verified, not assumed --
# see docs/otilib_discovery_report.md.

option(WARP_HYPAD_WITH_OTILIB "Build the OTIlib (GPL-3.0) sensitivity backend" OFF)

# The pin. OTIlib is a pinned external dependency, NOT vendored source: the
# checkout is git-ignored (see .gitignore, and docs/otilib_discovery_report.md
# for why: GPL-3.0 is accepted, and the checkout is kept external so a large
# upstream history stays out of this repository. Fetch it with
# scripts/bootstrap_otilib.sh or .ps1).
#
# The consequence of ignoring it is that the pin used to live only in a comment,
# so a silently different checkout would have produced a green gate against the
# wrong dependency. The revision is therefore verified at configure time.
set(WARP_HYPAD_OTILIB_PIN "de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54")

if(WARP_HYPAD_WITH_OTILIB)
  set(OTILIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/external/otilib")
  if(NOT EXISTS "${OTILIB_ROOT}/include/oti/static/onumm3n2.h")
    message(FATAL_ERROR
      "WARP_HYPAD_WITH_OTILIB=ON but OTIlib is missing at ${OTILIB_ROOT}.\n"
      "  git clone https://github.com/mauriaristi/otilib.git external/otilib\n"
      "  git -C external/otilib checkout ${WARP_HYPAD_OTILIB_PIN}")
  endif()

  # Verify the pin. A wrong revision is a FATAL_ERROR rather than a warning: a
  # sensitivity gate that passes against an unknown OTIlib revision is not a
  # gate. An absent git history is reported and tolerated, since a source
  # tarball is a legitimate way to supply the dependency.
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${OTILIB_ROOT}"
      OUTPUT_VARIABLE _oti_head
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE _oti_head_result)
    if(_oti_head_result EQUAL 0 AND _oti_head)
      if(NOT _oti_head STREQUAL WARP_HYPAD_OTILIB_PIN)
        message(FATAL_ERROR
          "OTIlib is at the wrong revision.\n"
          "  expected ${WARP_HYPAD_OTILIB_PIN}\n"
          "  found    ${_oti_head}\n"
          "  git -C external/otilib checkout ${WARP_HYPAD_OTILIB_PIN}")
      endif()
      message(STATUS "  otilib pin verified: ${WARP_HYPAD_OTILIB_PIN}")
    else()
      message(STATUS
        "  otilib pin NOT verified: no git history at ${OTILIB_ROOT}. "
        "Expected ${WARP_HYPAD_OTILIB_PIN}.")
    endif()
  endif()

  # The OTI C sources are third-party and are NOT held to this project's
  # warnings-as-errors policy; doing so would mean patching upstream, which
  # would fork the dependency.
  add_library(warp_hypad_otilib STATIC
    src/sensitivity/oti_onumm3n2_unity.c
    src/sensitivity/otilib_adapter.cpp
    "${OTILIB_ROOT}/src/c/real.c"
  )
  target_include_directories(warp_hypad_otilib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${OTILIB_ROOT}/include>
    $<BUILD_INTERFACE:${OTILIB_ROOT}/src/c>
  )
  target_link_libraries(warp_hypad_otilib PUBLIC warp_hypad_core)
  if(MSVC)
    target_compile_options(warp_hypad_otilib PRIVATE /w /fp:precise)
  else()
    target_compile_options(warp_hypad_otilib PRIVATE -w -ffp-contract=off)
  endif()

  add_executable(otilib_verify
    tests/otilib_verify.cpp
    tests/otilib_checks.cpp
    tests/otilib_fd_checks.cpp
  )
  target_link_libraries(otilib_verify PRIVATE warp_hypad_otilib)
  target_include_directories(otilib_verify PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  if(MSVC)
    target_compile_options(otilib_verify PRIVATE /W4 /permissive- /Zc:__cplusplus)
  else()
    target_compile_options(otilib_verify PRIVATE -Wall -Wextra)
  endif()

  set(WARP_HYPAD_OTI_RESULTS "${CMAKE_BINARY_DIR}/results/verification/otilib")
  add_test(NAME otilib.gate
           COMMAND otilib_verify --out "${WARP_HYPAD_OTI_RESULTS}")

  # Let A3 and A2T use the real HyPAD scalar for their sensitivity checks.
  target_link_libraries(milestone3_verify PRIVATE warp_hypad_otilib)
  target_compile_definitions(milestone3_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
  # The sweep exporter uses the OTI backend for sensitivities when it is
  # available, and centred finite differences on the same production metric when
  # it is not. Only the OTI COMPARISON cases skip in the disabled build; every
  # field is still exported either way.
  target_link_libraries(a_sweep_export PRIVATE warp_hypad_otilib)
  target_compile_definitions(a_sweep_export PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
  target_link_libraries(milestone2t_verify PRIVATE warp_hypad_otilib)
  target_compile_definitions(milestone2t_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
  if(TARGET b1_verify)
    target_link_libraries(b1_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b1_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
    target_link_libraries(b1c_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b1c_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
    target_link_libraries(b2_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b2_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
    target_link_libraries(b3_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b3_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
    target_link_libraries(b4_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b4_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
    target_link_libraries(b4c_verify PRIVATE warp_hypad_otilib)
    target_compile_definitions(b4c_verify PRIVATE WARP_HYPAD_HAVE_OTILIB=1)
  endif()

  message(STATUS "  otilib:   ENABLED (GPL-3.0 -- see docs/otilib_discovery_report.md)")
else()
  target_link_libraries(milestone3_verify PRIVATE warp_hypad_core)
  target_link_libraries(milestone2t_verify PRIVATE warp_hypad_core)
  message(STATUS "  otilib:   disabled (-DWARP_HYPAD_WITH_OTILIB=ON to enable)")
endif()

message(STATUS "warp_hypad ${PROJECT_VERSION}")
message(STATUS "  commit:   ${WARP_HYPAD_GIT_COMMIT} (dirty: ${WARP_HYPAD_GIT_DIRTY})")
message(STATUS "  compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  results:  ${WARP_HYPAD_RESULTS_DIR}")

# --- Program B figures ------------------------------------------------------
#
# Registered last, deliberately: it needs BOTH b1_verify and the Python3 lookup,
# and the Program B block is configured before the figures block. Placing it
# here rather than reordering keeps each section self-contained.

if(TARGET b1_verify AND Python3_FOUND)
  add_test(NAME programb1.figures
           COMMAND ${Python3_EXECUTABLE}
                   "${PROGRAM_B_ROOT}/python/plot_b1.py"
                   "${WARP_HYPAD_B1_RESULTS}")
  set_tests_properties(programb1.figures PROPERTIES DEPENDS "programb1.gate")

  add_test(NAME programb2.figures
           COMMAND ${Python3_EXECUTABLE}
                   "${PROGRAM_B_ROOT}/python/plot_b2.py"
                   "${WARP_HYPAD_B2_RESULTS}")
  set_tests_properties(programb2.figures PROPERTIES DEPENDS "programb2.gate")

  add_test(NAME programb3.figures
           COMMAND ${Python3_EXECUTABLE}
                   "${PROGRAM_B_ROOT}/python/plot_b3.py"
                   "${WARP_HYPAD_B3_RESULTS}")
  set_tests_properties(programb3.figures PROPERTIES DEPENDS "programb3.gate")

  add_test(NAME programb4.figures
           COMMAND ${Python3_EXECUTABLE}
                   "${PROGRAM_B_ROOT}/python/plot_b4.py"
                   "${WARP_HYPAD_B4_RESULTS}")
  set_tests_properties(programb4.figures PROPERTIES DEPENDS "programb4.gate")
endif()
