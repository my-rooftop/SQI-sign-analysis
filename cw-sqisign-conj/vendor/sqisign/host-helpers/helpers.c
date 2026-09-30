#include "quaternion_tests.h"
#include <quaternion_data.h>
#include <quaternion_constants.h>

// needed for the case where n_gamma is_even which is used in the tests
int
quat_represent_integer_even(quat_alg_elem_t *gamma, const ibz_t *n_gamma)
{
    // var dec
    int found;
    ibz_t cornacchia_target;
    ibz_t adjusted_n_gamma;
    ibz_t bound, sq_bound, temp;
    ibz_t test;
    ibz_vec_4_t coeffs; // coeffs = [x,y,z,t]
    quat_alg_elem_t quat_temp;

    // var init
    found = 0;
    ibz_init(&bound);
    ibz_init(&test);
    ibz_init(&temp);
    ibz_init(&sq_bound);
    ibz_vec_4_init(&coeffs);
    quat_alg_elem_init(&quat_temp);
    ibz_init(&adjusted_n_gamma);
    ibz_init(&cornacchia_target);

    // adjusting the norm of gamma (multiplying by 4 to find a solution in an order of odd level)
    ibz_mul_2exp(&adjusted_n_gamma, n_gamma, 2);
    // computation of the first bound = sqrt (adjust_n_gamma / p - 1)
    ibz_div(&sq_bound, &bound, &adjusted_n_gamma, &QUATALG_PINFTY.p);
    ibz_sub(&sq_bound, &sq_bound, &ibz_const_one);
    ibz_sqrt_floor(&bound, &sq_bound);
    ibz_set_bound(&bound, (ibz_get_bound(&adjusted_n_gamma) - ibz_bitsize(&QUATALG_PINFTY.p)) / 2 + 3);

    // the size of the search space is roughly n_gamma / p
    ibz_t counter;
    ibz_init(&counter);
    ibz_mul(&temp, &QUATALG_PINFTY.p, &QUATALG_PINFTY.p);
    ibz_sqrt_floor(&temp, &temp);
    ibz_div(&counter, &temp, &adjusted_n_gamma, &temp);
    ibz_set_bound(&counter, (int)((ibz_bitsize(&QUATALG_PINFTY.p) * 3) / 2));

    // entering the main loop
    while (!found && ibz_cmp(&counter, &ibz_const_zero) != 0) {
        // decreasing the counter
        ibz_sub(&counter, &counter, &ibz_const_one);
        ibz_set_bound(&counter, (int)((ibz_bitsize(&QUATALG_PINFTY.p) * 3) / 2));

        // we start by sampling the first coordinate
        if (!ibz_rand_interval(&coeffs.v[2], &ibz_const_one, &bound))
            break;

        // then, we sample the second coordinate
        // computing the second bound in temp as sqrt( (adjust_n_gamma - p*coeffs.v[2]²)/qp )
        ibz_mul(&cornacchia_target, &coeffs.v[2], &coeffs.v[2]);
        ibz_mul(&temp, &cornacchia_target, &QUATALG_PINFTY.p);
        ibz_sub(&temp, &adjusted_n_gamma, &temp);
        ibz_div(&temp, &sq_bound, &temp, &QUATALG_PINFTY.p);
        ibz_set_bound(&temp, ibz_get_bound(&adjusted_n_gamma) - ibz_bitsize(&QUATALG_PINFTY.p) + 2);
        ibz_sqrt_floor(&temp, &temp);

        if (ibz_cmp(&temp, &ibz_const_zero) <= 0) {
            continue;
        }
        // sampling the second value
        if (!ibz_rand_interval(&coeffs.v[3], &ibz_const_one, &temp))
            break;

        // compute cornacchia_target = n_gamma - p * (z² + t²)
        ibz_mul(&temp, &coeffs.v[3], &coeffs.v[3]);
        ibz_add(&cornacchia_target, &cornacchia_target, &temp);
        ibz_mul(&cornacchia_target, &cornacchia_target, &QUATALG_PINFTY.p);
        ibz_sub(&cornacchia_target, &adjusted_n_gamma, &cornacchia_target);
        ibz_set_bound(&cornacchia_target, ibz_get_bound(n_gamma) + 2);
        assert(ibz_cmp(&cornacchia_target, &ibz_const_zero) > 0);

        // applying cornacchia
        if (ibz_probab_prime(&cornacchia_target, QUAT_primality_num_iter))
            found = ibz_cornacchia_prime(&(coeffs.v[0]), &(coeffs.v[1]), &cornacchia_target);
        else
            found = 0;

        if (found) {
            for (int i = 0; i < 2; i++)
                ibz_set_bound(&coeffs.v[i], (ibz_get_bound(&cornacchia_target) / 2) + 2);

#ifndef NDEBUG
            ibz_set(&temp, 1, 2);
            ibz_mul(&temp, &temp, &(coeffs.v[1]));
            ibz_mul(&temp, &temp, &(coeffs.v[1]));
            ibz_mul(&test, &(coeffs.v[0]), &(coeffs.v[0]));
            ibz_add(&temp, &temp, &test);
            assert(0 == ibz_cmp(&temp, &cornacchia_target));

            ibz_mul(&cornacchia_target, &(coeffs.v[3]), &(coeffs.v[3]));
            ibz_mul(&cornacchia_target, &cornacchia_target, &QUATALG_PINFTY.p);
            ibz_mul(&temp, &(coeffs.v[1]), &(coeffs.v[1]));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs.v[0]), &coeffs.v[0]);
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs.v[2]), &coeffs.v[2]);
            ibz_mul(&temp, &temp, &QUATALG_PINFTY.p);
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            assert(0 == ibz_cmp(&cornacchia_target, &adjusted_n_gamma));
#endif
            // translate x,y,z,t into the quaternion element gamma, assuming 1,i,j,ij basis
            ibz_vec_4_copy(&gamma->coord, &coeffs);
            ibz_set(&gamma->denom, 1, 2);
#ifndef NDEBUG
            quat_alg_norm(&temp, &(coeffs.v[0]), gamma, &QUATALG_PINFTY);
            assert(ibz_is_one(&(coeffs.v[0])));
            assert(0 == ibz_cmp(&temp, &adjusted_n_gamma));
            assert(quat_lattice_contains(NULL, &MAXORD_O0, gamma));
#endif
            // making gamma primitive
            // coeffs contains the coefficients of primitivized gamma in the basis of order
            quat_alg_make_primitive(&coeffs, &temp, gamma, &MAXORD_O0);

            found = (ibz_cmp(&temp, &ibz_const_two) == 0);
        }
    }

    if (found) {
        // new gamma
        ibz_mat_4x4_eval(&coeffs, &MAXORD_O0.basis, &coeffs);
        ibz_copy(&gamma->coord.v[0], &coeffs.v[0]);
        ibz_copy(&gamma->coord.v[1], &coeffs.v[1]);
        ibz_copy(&gamma->coord.v[2], &coeffs.v[2]);
        ibz_copy(&gamma->coord.v[3], &coeffs.v[3]);
        ibz_copy(&gamma->denom, &MAXORD_O0.denom);
    }
    return (found);
}

int
quat_lattice_is_ideal(const quat_lattice_t *lat)
{
    int res = 1;
    quat_alg_elem_t elem, a, b;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&elem);
    ibz_copy(&a.denom, &lat->denom);
    ibz_copy(&b.denom, &lat->denom);
    for (int i = 0; i < 4; i++) {
        ibz_vec_4_copy_ibz(
            &a.coord, &lat->basis.m[0][i], &lat->basis.m[1][i], &lat->basis.m[2][i], &lat->basis.m[3][i]);
        for (int j = 0; j < 4; j++) {
            ibz_vec_4_copy_ibz(
                &b.coord, &lat->basis.m[0][j], &lat->basis.m[1][j], &lat->basis.m[2][j], &lat->basis.m[3][j]);
            quat_alg_mul(&elem, &a, &b, &QUATALG_PINFTY);
            res = res && quat_lattice_contains(NULL, lat, &elem);
        }
    }
    return (res);
}

int
quat_lattice_is_O0_ideal(const quat_lattice_t *lat)
{
    int res = quat_lattice_is_ideal(lat);
    quat_alg_elem_t elem, a, b;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&elem);
    ibz_copy(&a.denom, &MAXORD_O0.denom);
    ibz_copy(&b.denom, &lat->denom);
    for (int i = 0; i < 4; i++) {
        ibz_vec_4_copy_ibz(&a.coord,
                           &MAXORD_O0.basis.m[0][i],
                           &MAXORD_O0.basis.m[1][i],
                           &MAXORD_O0.basis.m[2][i],
                           &MAXORD_O0.basis.m[3][i]);
        for (int j = 0; j < 4; j++) {
            ibz_vec_4_copy_ibz(
                &b.coord, &lat->basis.m[0][j], &lat->basis.m[1][j], &lat->basis.m[2][j], &lat->basis.m[3][j]);
            quat_alg_mul(&elem, &a, &b, &QUATALG_PINFTY);
            res = res && quat_lattice_contains(NULL, lat, &elem);
        }
    }
    return (res);
}

int
quat_lattice_is_maximal_order(const quat_lattice_t *lat)
{
    ibz_t norm;
    ibz_init(&norm);
    int res = quat_lattice_is_ideal(lat);
    if (res) {
        quat_lattice_norm(&norm, lat);
        res = res && ibz_is_one(&norm);
    }
    return (res);
}

void
quat_lattice_norm(ibz_t *norm, const quat_lattice_t *ideal)
{
    ibz_t n;
    ibz_init(&n);
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &n, &ideal->basis);
    ibz_add(&n, &n, &n);
    ibz_add(norm, &n, &n);
    ibz_mul(&n, &ideal->denom, &ideal->denom);
    ibz_mul(&n, &n, &n);
    ibz_div(norm, &n, norm, &n);
    ibz_abs(norm, norm);
    assert(ibz_is_zero(&n));
    ibz_sqrt_floor(norm, norm);
}

void
quat_lattice_from_ideal_split(quat_lattice_t *lat, const quat_ideal_t *ideal, const quat_alg_elem_t *split)
{
    quat_alg_elem_t elem;
    quat_lattice_t start;
    quat_lattice_init(&start);
    quat_alg_elem_init(&elem);
    quat_to_lattice(&start, ideal);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&lat->denom, &start.denom, &split->denom);
        ibz_vec_4_copy_ibz(
            &elem.coord, &start.basis.m[0][i], &start.basis.m[1][i], &start.basis.m[2][i], &start.basis.m[3][i]);
        quat_alg_coord_mul(&elem.coord, &elem.coord, &split->coord, &QUATALG_PINFTY);

        ibz_copy(&lat->basis.m[0][i], &elem.coord.v[0]);
        ibz_copy(&lat->basis.m[1][i], &elem.coord.v[1]);
        ibz_copy(&lat->basis.m[2][i], &elem.coord.v[2]);
        ibz_copy(&lat->basis.m[3][i], &elem.coord.v[3]);
    }
    quat_lattice_reduce_denom(lat, lat);
}

int
quat_ideal_shape(const quat_ideal_t *ideal)
{
    int res = 0;
    res = res || !(ibz_cmp(&ideal->x, &ideal->norm) < 0);
    res = res || !(ibz_cmp(&ideal->y, &ideal->norm) < 0);
    res = res || !(ibz_cmp(&ideal->x, &ibz_const_zero) >= 0);
    res = res || !(ibz_cmp(&ideal->y, &ibz_const_zero) >= 0);
    return (!res);
}

int
quat_ideal_equal(const quat_ideal_t *a, const quat_ideal_t *b)
{
    int res = 1;
    res = res && (ibz_cmp(&a->x, &b->x) == 0);
    res = res && (ibz_cmp(&a->y, &b->y) == 0);
    res = res && (ibz_cmp(&a->norm, &b->norm) == 0);
    return (res);
}

int
quat_ideal_is_ideal(const quat_ideal_t *ideal)
{
    quat_lattice_t lat;
    quat_lattice_init(&lat);
    quat_to_lattice(&lat, ideal);
    int res = quat_lattice_is_O0_ideal(&lat);
    return (res);
}

int
quat_ideals_equivalence(const quat_ideal_t *a, const quat_ideal_t *b)
{
    int res = 1;
    quat_ideal_t a_red, b_red;
    quat_ideal_init(&a_red);
    quat_ideal_init(&b_red);
    quat_ideal_shortest_equivalent(NULL, NULL, &b_red, b);
    quat_ideal_shortest_equivalent(NULL, NULL, &a_red, a);
    res = quat_ideal_equal(&b_red, &a_red);
    return (res);
}
