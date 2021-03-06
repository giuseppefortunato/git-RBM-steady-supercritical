#project(utopiae-orbit-dependencies
#        LANGUAGES CXX
#        VERSION 1.0.0)

# Include external libraries
include(ExternalProject)

#set(EIGEN_INCLUDE_DIR /usr/include/eigen3/)
set(EIGEN_INCLUDE_DIR ${PROJECT_PATH}/include/Eigen)
include_directories (${EIGEN_INCLUDE_DIR})


# When gsl is required
if (ENABLE_GSL)
    find_package(GSL REQUIRED)
    list(APPEND MANDATORY_LIBRARIES "${GSL_LIBRARIES}")
    include_directories(${GSL_INCLUDE_DIRS})
endif()


# Pagmo is optional
if (ENABLE_PAGMO)

    # Add compiler flag
    add_definitions(-D__USE_PAGMO)

    include(cmake/external/pagmo.cmake)
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_STATIC_LIBRARY}")
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_LIBRARY}")

    include_directories (${PAGMO_INCLUDE_DIR})

endif()


