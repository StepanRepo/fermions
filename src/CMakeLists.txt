add_library(Precision precision.f95)

# Define mathematical functions
# needed to compute integrals with gaussian quadratures
add_library(SLE sle_solve.f95)
add_dependencies(SLE Precision)
target_link_libraries(SLE PUBLIC Precision)

add_library(Polynomials polynomials.f95)
add_dependencies(Polynomials Precision)
target_link_libraries(Polynomials PUBLIC Precision)

add_library(GaussInt gauss_int.f95)
add_dependencies(GaussInt SLE Precision Polynomials)
target_link_libraries(GaussInt PUBLIC SLE Precision Polynomials)


add_library(Distributions distributions.f95)
add_dependencies(Distributions GaussInt Precision Config)
target_link_libraries(Distributions PUBLIC GaussInt Precision Config)


# Define configuration
add_library(Config config.f95)
add_dependencies(Config Precision)


# Make final executable
add_executable(fermions main.f95)

target_link_libraries(fermions PUBLIC 
	Precision GaussInt Config Polynomials SLE Distributions
	LAPACK::LAPACK)


# Move executable to thr working directiory
install(TARGETS fermions DESTINATION ${CMAKE_BIN_DIR})
