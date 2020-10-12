#CMake file for Adding executables

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    if (BUILD_STATIC)
        target_link_libraries(${arg1} MODES_lib_static ${MANDATORY_LIBRARIES})
        #target_link_libraries(${arg1} PUBLIC "-static")
    else ()
        target_link_libraries(${arg1} MODES ${MANDATORY_LIBRARIES})
    endif()
endfunction()

#ADD_RBM_TESTCASE(AdaptiveMODES_Offline)
#ADD_RBM_TESTCASE(AdaptiveMODES_Online)
#ADD_RBM_TESTCASE(AdaptiveMODES_DTR_Offline)
#ADD_RBM_TESTCASE(AdaptiveMODES_DTR_Online)
ADD_RBM_TESTCASE(SingleMODES)
#ADD_RBM_TESTCASE(test_rbf)
#ADD_RBM_TESTCASE(POD_Rec)
#ADD_RBM_TESTCASE(Stesttest)


if(ENABLE_PAGMO)
#    ADD_RBM_TESTCASE(Adaptive_Sampling)
#    ADD_RBM_TESTCASE(AdaptSample_ResEval)
endif()

