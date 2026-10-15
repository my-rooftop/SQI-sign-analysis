#ifndef CW_SQISIGN_RESPONSE_CORE_H
#define CW_SQISIGN_RESPONSE_CORE_H

#include <stdint.h>

#include <quaternion.h>

typedef struct {
    uint32_t parse_cycles;
    uint32_t core_cycles;
    uint32_t serialize_cycles;
    uint32_t handler_cycles;
    uint32_t inter_cycles;
    uint32_t conjugate_cycles;
    uint32_t lattice_intersect_cycles;
    uint32_t content_mul_cycles;
    uint32_t sample_total_cycles;
    uint32_t sample_setup_cycles;
    uint32_t sample_rejection_cycles;
    uint32_t sample_finalize_cycles;
    uint32_t normalize_total_cycles;
    uint32_t normalize_gcd_cycles;
    uint32_t normalize_coord0_div_cycles;
    uint32_t normalize_coord_rest_div_cycles;
    uint32_t normalize_sign_cycles;
    uint32_t rejection_iters;
} response_profile_t;

extern volatile response_profile_t response_profile;

void response_profile_reset(void);
void response_cycles_init(void);
uint32_t response_cycles_read(void);

void sqisign_response_quat_element_core(quat_alg_elem_t *resp_quat,
                                        ibz_t *lattice_content,
                                        const quat_left_ideal_t *secret_ideal,
                                        const quat_left_ideal_t *challenge_ideal,
                                        const quat_left_ideal_t *commit_ideal);

#endif
