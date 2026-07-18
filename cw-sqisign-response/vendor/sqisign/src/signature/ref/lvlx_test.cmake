add_executable(sqisign_test_signature_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_signature.c)
target_link_libraries(sqisign_test_signature_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_test)
target_include_directories(sqisign_test_signature_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_executable(sqisign_test_threadsafety_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_threadsafety.c)
target_link_libraries(sqisign_test_threadsafety_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_test pthread)
target_include_directories(sqisign_test_threadsafety_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_test(sqisign_test_signature_${SVARIANT_LOWER} sqisign_test_signature_${SVARIANT_LOWER} 3)
add_test(sqisign_test_threadsafety_${SVARIANT_LOWER} sqisign_test_threadsafety_${SVARIANT_LOWER} 3)

add_custom_command(
  TARGET sqisign_test_signature_${SVARIANT_LOWER}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy $<TARGET_FILE:sqisign_test_signature_${SVARIANT_LOWER}> ${CMAKE_BINARY_DIR}/test/sqisign_${SVARIANT_LOWER}
)

add_executable(sqisign_bench_signature_${SVARIANT_LOWER} ${LVLX_DIR}/test/bench_signature.c)
target_link_libraries(sqisign_bench_signature_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_sys)
target_include_directories(sqisign_bench_signature_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_executable(sqisign_check_dumped_key_${SVARIANT_LOWER} ${LVLX_DIR}/test/check_dumped_key.c)
target_link_libraries(sqisign_check_dumped_key_${SVARIANT_LOWER} sqisign_${SVARIANT_LOWER}_test_nistapi sqisign_common_sys)
target_include_directories(sqisign_check_dumped_key_${SVARIANT_LOWER} PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_executable(sqisign_check_dumped_key_lcg_${SVARIANT_LOWER} ${LVLX_DIR}/test/check_dumped_key_lcg.c)
target_compile_definitions(sqisign_check_dumped_key_lcg_${SVARIANT_LOWER} PRIVATE SQISIGN_CW_PROGRESS NDEBUG)
target_compile_definitions(${LIB_SIGNATURE_${SVARIANT_UPPER}} PRIVATE SQISIGN_CW_PROGRESS NDEBUG)
target_compile_definitions(${LIB_ID2ISO_${SVARIANT_UPPER}} PRIVATE SQISIGN_CW_PROGRESS NDEBUG)
target_compile_definitions(${LIB_HD_${SVARIANT_UPPER}} PRIVATE SQISIGN_CW_PROGRESS NDEBUG)
target_compile_definitions(${LIB_QUATERNION} PRIVATE SQISIGN_CW_PROGRESS NDEBUG)
target_link_libraries(sqisign_check_dumped_key_lcg_${SVARIANT_LOWER} sqisign_${SVARIANT_LOWER}_test_nistapi)
target_include_directories(sqisign_check_dumped_key_lcg_${SVARIANT_LOWER} PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_executable(sqisign_ram_budget_signature_${SVARIANT_LOWER} ${LVLX_DIR}/test/bench_signature.c ${LVLX_DIR}/test/ram_budget_probe.c)
target_compile_definitions(sqisign_ram_budget_signature_${SVARIANT_LOWER} PRIVATE SQISIGN_CW_PROGRESS)
target_compile_definitions(${LIB_SIGNATURE_${SVARIANT_UPPER}} PRIVATE SQISIGN_CW_PROGRESS)
target_compile_definitions(${LIB_ID2ISO_${SVARIANT_UPPER}} PRIVATE SQISIGN_CW_PROGRESS)
target_link_libraries(sqisign_ram_budget_signature_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_sys)
target_include_directories(sqisign_ram_budget_signature_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})
target_link_options(sqisign_ram_budget_signature_${SVARIANT_LOWER} PRIVATE -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=free)

set(BM_BINS ${BM_BINS} sqisign_bench_signature_${SVARIANT_LOWER} CACHE INTERNAL "List of benchmark executables")
