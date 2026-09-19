
add_executable(test_wgt_storage test_wgt_storage.c)
# GMP arrives through edge_weight_storage, which links GMP::GMP publicly for
# qisq2; naming it again as a bare -lgmp would miss a non-default prefix.
target_link_libraries(test_wgt_storage PRIVATE edge_weight_storage m)

#add_executable(test_mpfr_map test_mpfr_map.c)
#target_link_libraries(test_mpfr_map edge_weight_storage)

add_test(NAME test_wgt_storage COMMAND test_wgt_storage)
#add_test(test_mpfr_map test_mpfr_map)