#include <assert.h>
#include <bench.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <encoded_sizes.h>
#include <id2iso.h>
#include <quaternion_data.h>
#include <rng.h>
#include <signature.h>
#include <torsion_constants.h>
#include <tools.h>

#include "intbig_internal.h"

static bool commit(ec_curve_t *E_com, ec_basis_t *basis_even_com, quat_left_ideal_t *lideal_com)
{
    bool found = quat_sampling_random_ideal_O0_given_norm(
        lideal_com, &COM_DEGREE, 1, &QUAT_represent_integer_params, NULL);
    found = found && quat_lideal_prime_norm_reduced_equivalent(
                         lideal_com, &QUATALG_PINFTY, QUAT_primality_num_iter, QUAT_equiv_bound_coeff);
    found = found && dim2id2iso_arbitrary_isogeny_evaluation(basis_even_com, E_com, lideal_com);
    return found;
}

static void compute_challenge_ideal_signature(quat_left_ideal_t *lideal_chall_two,
                                              const signature_t *sig,
                                              const secret_key_t *sk)
{
    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);

    ibz_set(&vec[0], 1);
    ibz_copy_digit_array(&vec[1], sig->chall_coeff);
    ibz_mat_2x2_eval(&vec, &(sk->mat_BAcan_to_BA0_two), &vec);
    id2iso_kernel_dlogs_to_ideal_even(lideal_chall_two, &vec, TORSION_EVEN_POWER);

    ibz_vec_2_finalize(&vec);
}

static void die(const char *msg)
{
    fprintf(stderr, "error: %s\n", msg);
    exit(1);
}

static char *read_all(FILE *f)
{
    size_t cap = 4096;
    size_t len = 0;
    char *buf = malloc(cap);
    if (!buf) {
        die("malloc failed");
    }
    for (;;) {
        if (len + 2048 + 1 > cap) {
            cap *= 2;
            char *tmp = realloc(buf, cap);
            if (!tmp) {
                free(buf);
                die("realloc failed");
            }
            buf = tmp;
        }
        size_t n = fread(buf + len, 1, 2048, f);
        len += n;
        if (n < 2048) {
            if (ferror(f)) {
                free(buf);
                die("read failed");
            }
            break;
        }
    }
    buf[len] = '\0';
    return buf;
}

static char *next_token(char **cursor)
{
    char *s = *cursor;
    while (*s == '\n' || *s == '\r' || *s == ' ' || *s == '\t') {
        s++;
    }
    if (*s == '\0') {
        *cursor = s;
        return NULL;
    }
    char *tok = s;
    while (*s != '\0' && *s != '\n' && *s != '\r' && *s != ' ' && *s != '\t') {
        s++;
    }
    if (*s != '\0') {
        *s++ = '\0';
    }
    *cursor = s;
    return tok;
}

static void parse_ibz_or_die(ibz_t *x, char **cursor)
{
    char *tok = next_token(cursor);
    if (!tok || !ibz_set_from_str(x, tok, 16)) {
        die("bad integer input");
    }
}

static void print_ibz_line(const ibz_t *x)
{
    char tmp[2048];
    if (!ibz_convert_to_str(x, tmp, 16)) {
        die("integer conversion failed");
    }
    puts(tmp);
}

static void print_ibz_inline(const ibz_t *x, char sep)
{
    char tmp[2048];
    if (!ibz_convert_to_str(x, tmp, 16)) {
        die("integer conversion failed");
    }
    fputs(tmp, stdout);
    putchar(sep);
}

static void print_ideal(const quat_left_ideal_t *ideal)
{
    print_ibz_line(&ideal->norm);
    print_ibz_line(&ideal->lattice.denom);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            print_ibz_line(&ideal->lattice.basis[i][j]);
        }
    }
}

static void sample_response_host(quat_alg_elem_t *x,
                                 const quat_lattice_t *lattice,
                                 const ibz_t *lattice_content)
{
    ibz_t bound;
    ibz_init(&bound);
    ibz_pow(&bound, &ibz_const_two, SQIsign_response_length);
    ibz_sub(&bound, &bound, &ibz_const_one);
    ibz_mul(&bound, &bound, lattice_content);

    int ok = quat_lattice_sample_from_ball(x, lattice, &QUATALG_PINFTY, &bound);
    if (!ok) {
        die("response sampling failed");
    }

    ibz_finalize(&bound);
}

static void compute_response_quat_element_host(quat_alg_elem_t *resp_quat,
                                               ibz_t *lattice_content,
                                               const secret_key_t *sk,
                                               const quat_left_ideal_t *challenge_ideal,
                                               const quat_left_ideal_t *commit_ideal)
{
    quat_left_ideal_t challenge_secret;
    quat_lattice_t commit_conj;
    quat_lattice_t hom_lattice;

    quat_left_ideal_init(&challenge_secret);
    quat_lattice_init(&commit_conj);
    quat_lattice_init(&hom_lattice);

    quat_lideal_inter(&challenge_secret, challenge_ideal, &sk->secret_ideal, &QUATALG_PINFTY);
    quat_lattice_conjugate_without_hnf(&commit_conj, &commit_ideal->lattice);
    quat_lattice_intersect(&hom_lattice, &challenge_secret.lattice, &commit_conj);
    ibz_mul(lattice_content, &challenge_secret.norm, &commit_ideal->norm);
    sample_response_host(resp_quat, &hom_lattice, lattice_content);

    quat_lattice_finalize(&hom_lattice);
    quat_lattice_finalize(&commit_conj);
    quat_left_ideal_finalize(&challenge_secret);
}

static void compute_response_quat_element_from_ideals(quat_alg_elem_t *resp_quat,
                                                      ibz_t *lattice_content,
                                                      const quat_left_ideal_t *secret_ideal,
                                                      const quat_left_ideal_t *challenge_ideal,
                                                      const quat_left_ideal_t *commit_ideal)
{
    quat_left_ideal_t challenge_secret;
    quat_lattice_t commit_conj;
    quat_lattice_t hom_lattice;

    quat_left_ideal_init(&challenge_secret);
    quat_lattice_init(&commit_conj);
    quat_lattice_init(&hom_lattice);

    quat_lideal_inter(&challenge_secret, challenge_ideal, secret_ideal, &QUATALG_PINFTY);
    quat_lattice_conjugate_without_hnf(&commit_conj, &commit_ideal->lattice);
    quat_lattice_intersect(&hom_lattice, &challenge_secret.lattice, &commit_conj);
    ibz_mul(lattice_content, &challenge_secret.norm, &commit_ideal->norm);
    sample_response_host(resp_quat, &hom_lattice, lattice_content);

    quat_lattice_finalize(&hom_lattice);
    quat_lattice_finalize(&commit_conj);
    quat_left_ideal_finalize(&challenge_secret);
}

static void parse_ideal_or_die(quat_left_ideal_t *ideal, char **cursor)
{
    ideal->parent_order = &MAXORD_O0;
    parse_ibz_or_die(&ideal->norm, cursor);
    parse_ibz_or_die(&ideal->lattice.denom, cursor);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            parse_ibz_or_die(&ideal->lattice.basis[i][j], cursor);
        }
    }
}

static void seed_rng(void)
{
    unsigned char seed[48];
    for (size_t i = 0; i < sizeof(seed); i++) {
        seed[i] = (unsigned char)i;
    }
    randombytes_init(seed, NULL, 256);
}

static int cmd_gen(void)
{
    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    ec_curve_t Ecom;
    ec_basis_t basis_even_com;
    quat_left_ideal_t lideal_commit;
    quat_left_ideal_t lideal_chall_two;
    unsigned char msg[32] = { 0 };

    memset(&pk, 0, sizeof(pk));
    memset(&sig, 0, sizeof(sig));
    secret_key_init(&sk);
    ec_curve_init(&Ecom);
    quat_left_ideal_init(&lideal_commit);
    quat_left_ideal_init(&lideal_chall_two);

    seed_rng();
    if (!protocols_keygen(&pk, &sk)) {
        die("keygen failed");
    }
    if (!commit(&Ecom, &basis_even_com, &lideal_commit)) {
        die("commit failed");
    }
    hash_to_challenge(&sig.chall_coeff, &pk, &Ecom, msg, sizeof(msg));
    compute_challenge_ideal_signature(&lideal_chall_two, &sig, &sk);

    print_ideal(&sk.secret_ideal);
    print_ideal(&lideal_chall_two);
    print_ideal(&lideal_commit);

    quat_left_ideal_finalize(&lideal_chall_two);
    quat_left_ideal_finalize(&lideal_commit);
    secret_key_finalize(&sk);
    return 0;
}

static unsigned ibz_bitlen_via_string(const ibz_t *x)
{
    char tmp[2048];
    if (!ibz_convert_to_str(x, tmp, 16)) {
        die("integer conversion failed");
    }
    char *s = tmp;
    if (*s == '-') {
        s++;
    }
    while (*s == '0' && s[1] != '\0') {
        s++;
    }
    if (strcmp(s, "0") == 0) {
        return 0;
    }
    unsigned hex_digits = (unsigned)strlen(s);
    unsigned top = 0;
    char c = s[0];
    if (c >= '0' && c <= '9') {
        top = (unsigned)(c - '0');
    } else if (c >= 'a' && c <= 'f') {
        top = (unsigned)(c - 'a' + 10);
    } else if (c >= 'A' && c <= 'F') {
        top = (unsigned)(c - 'A' + 10);
    }
    unsigned top_bits = 0;
    while (top) {
        top_bits++;
        top >>= 1;
    }
    return 4 * (hex_digits - 1) + top_bits;
}

static uint32_t ibz_low32_via_string(const ibz_t *x)
{
    char tmp[2048];
    if (!ibz_convert_to_str(x, tmp, 16)) {
        die("integer conversion failed");
    }
    char *s = tmp;
    if (*s == '-') {
        s++;
    }
    size_t len = strlen(s);
    uint32_t out = 0;
    size_t start = len > 8 ? len - 8 : 0;
    for (size_t i = start; i < len; i++) {
        char c = s[i];
        unsigned v = 0;
        if (c >= '0' && c <= '9') {
            v = (unsigned)(c - '0');
        } else if (c >= 'a' && c <= 'f') {
            v = (unsigned)(c - 'a' + 10);
        } else if (c >= 'A' && c <= 'F') {
            v = (unsigned)(c - 'A' + 10);
        }
        out = (out << 4) | v;
    }
    return out;
}

static void update_range(unsigned *min_bits, unsigned *max_bits, const ibz_t *x)
{
    unsigned bits = ibz_bitlen_via_string(x);
    if (bits < *min_bits) {
        *min_bits = bits;
    }
    if (bits > *max_bits) {
        *max_bits = bits;
    }
}

static void print_value_summary(const char *name, const ibz_t *x)
{
    printf("%s_sign=%d %s_bits=%u %s_limb_count32=%u %s_low32=%08x ",
           name,
           ibz_cmp(x, &ibz_const_zero) < 0,
           name,
           ibz_bitlen_via_string(x),
           name,
           (ibz_bitlen_via_string(x) + 31) / 32,
           name,
           ibz_low32_via_string(x));
}

static int cmd_sample(unsigned long count)
{
    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    ec_curve_t Ecom;
    ec_basis_t basis_even_com;
    quat_left_ideal_t lideal_commit;
    quat_left_ideal_t lideal_chall_two;
    quat_alg_elem_t resp;
    ibz_t lattice_content;
    unsigned min_bits[6] = { 9999, 9999, 9999, 9999, 9999, 9999 };
    unsigned max_bits[6] = { 0, 0, 0, 0, 0, 0 };
    unsigned char msg[32] = { 0 };

    seed_rng();
    printf("idx ");
    printf("lattice_content_sign lattice_content_bits lattice_content_limb_count32 lattice_content_low32 ");
    printf("denom_sign denom_bits denom_limb_count32 denom_low32 ");
    for (int i = 0; i < 4; i++) {
        printf("coord%d_sign coord%d_bits coord%d_limb_count32 coord%d_low32 ", i, i, i, i);
    }
    printf("\n");

    for (unsigned long sample = 0; sample < count; sample++) {
        memset(&pk, 0, sizeof(pk));
        memset(&sig, 0, sizeof(sig));
        secret_key_init(&sk);
        ec_curve_init(&Ecom);
        quat_left_ideal_init(&lideal_commit);
        quat_left_ideal_init(&lideal_chall_two);
        quat_alg_elem_init(&resp);
        ibz_init(&lattice_content);

        if (!protocols_keygen(&pk, &sk)) {
            die("keygen failed");
        }
        if (!commit(&Ecom, &basis_even_com, &lideal_commit)) {
            die("commit failed");
        }
        hash_to_challenge(&sig.chall_coeff, &pk, &Ecom, msg, sizeof(msg));
        compute_challenge_ideal_signature(&lideal_chall_two, &sig, &sk);
        compute_response_quat_element_host(&resp, &lattice_content, &sk, &lideal_chall_two, &lideal_commit);

        update_range(&min_bits[0], &max_bits[0], &lattice_content);
        update_range(&min_bits[1], &max_bits[1], &resp.denom);
        for (int i = 0; i < 4; i++) {
            update_range(&min_bits[2 + i], &max_bits[2 + i], &resp.coord[i]);
        }

        printf("%lu ", sample);
        print_value_summary("lattice_content", &lattice_content);
        print_value_summary("denom", &resp.denom);
        for (int i = 0; i < 4; i++) {
            char name[16];
            snprintf(name, sizeof(name), "coord%d", i);
            print_value_summary(name, &resp.coord[i]);
        }
        printf("\n");

        ibz_finalize(&lattice_content);
        quat_alg_elem_finalize(&resp);
        quat_left_ideal_finalize(&lideal_chall_two);
        quat_left_ideal_finalize(&lideal_commit);
        secret_key_finalize(&sk);
    }

    const char *names[6] = { "lattice_content", "denom", "coord0", "coord1", "coord2", "coord3" };
    fprintf(stderr, "summary_bits");
    for (int i = 0; i < 6; i++) {
        fprintf(stderr, " %s=%u..%u", names[i], min_bits[i], max_bits[i]);
    }
    fprintf(stderr, "\n");
    return 0;
}

/* Compute N real responses and print each on one line as space-separated hex:
   "content denom coord0 coord1 coord2 coord3".  This is exactly the token order
   the board's cmd 'k' (response_set_known) parses, so the capture host can inject
   genuine response values -- real signs, real denom, real limb structure -- while
   still measuring only the fast serialization path (no ~12 s core on the board). */
static int cmd_dump_responses(unsigned long count)
{
    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    ec_curve_t Ecom;
    ec_basis_t basis_even_com;
    quat_left_ideal_t lideal_commit;
    quat_left_ideal_t lideal_chall_two;
    quat_alg_elem_t resp;
    ibz_t lattice_content;
    unsigned char msg[32] = { 0 };

    seed_rng();
    for (unsigned long sample = 0; sample < count; sample++) {
        memset(&pk, 0, sizeof(pk));
        memset(&sig, 0, sizeof(sig));
        secret_key_init(&sk);
        ec_curve_init(&Ecom);
        quat_left_ideal_init(&lideal_commit);
        quat_left_ideal_init(&lideal_chall_two);
        quat_alg_elem_init(&resp);
        ibz_init(&lattice_content);

        if (!protocols_keygen(&pk, &sk)) {
            die("keygen failed");
        }
        if (!commit(&Ecom, &basis_even_com, &lideal_commit)) {
            die("commit failed");
        }
        hash_to_challenge(&sig.chall_coeff, &pk, &Ecom, msg, sizeof(msg));
        compute_challenge_ideal_signature(&lideal_chall_two, &sig, &sk);
        compute_response_quat_element_host(&resp, &lattice_content, &sk, &lideal_chall_two, &lideal_commit);

        print_ibz_inline(&lattice_content, ' ');
        print_ibz_inline(&resp.denom, ' ');
        print_ibz_inline(&resp.coord[0], ' ');
        print_ibz_inline(&resp.coord[1], ' ');
        print_ibz_inline(&resp.coord[2], ' ');
        print_ibz_inline(&resp.coord[3], '\n');

        ibz_finalize(&lattice_content);
        quat_alg_elem_finalize(&resp);
        quat_left_ideal_finalize(&lideal_chall_two);
        quat_left_ideal_finalize(&lideal_commit);
        secret_key_finalize(&sk);

        if ((sample + 1) % 1000 == 0) {
            fprintf(stderr, "dumped %lu/%lu\n", sample + 1, count);
        }
    }
    return 0;
}

static void parse_response_or_die(quat_alg_elem_t *resp, ibz_t *lattice_content, char *text)
{
    char *cursor = text;
    char *status = next_token(&cursor);
    if (!status || strcmp(status, "1") != 0) {
        die("board response status is not 1");
    }
    parse_ibz_or_die(lattice_content, &cursor);
    parse_ibz_or_die(&resp->denom, &cursor);
    for (int i = 0; i < 4; i++) {
        parse_ibz_or_die(&resp->coord[i], &cursor);
    }
}

static int cmd_verify(const char *vector_path, const char *response_path)
{
    FILE *vf = fopen(vector_path, "rb");
    FILE *rf = fopen(response_path, "rb");
    if (!vf || !rf) {
        die("open input files failed");
    }
    char *vtext = read_all(vf);
    char *rtext = read_all(rf);
    fclose(vf);
    fclose(rf);

    char *cursor = vtext;
    quat_left_ideal_t secret_ideal;
    quat_left_ideal_t challenge_ideal;
    quat_left_ideal_t commit_ideal;
    quat_left_ideal_t challenge_secret;
    quat_lattice_t commit_conj;
    quat_lattice_t hom_lattice;
    quat_alg_elem_t resp;
    ibz_t lattice_content;
    ibz_t expected_content;
    ibz_t bound;
    ibz_t norm_num;
    ibz_t norm_den;
    ibz_t norm_bound;

    quat_left_ideal_init(&secret_ideal);
    quat_left_ideal_init(&challenge_ideal);
    quat_left_ideal_init(&commit_ideal);
    quat_left_ideal_init(&challenge_secret);
    quat_lattice_init(&commit_conj);
    quat_lattice_init(&hom_lattice);
    quat_alg_elem_init(&resp);
    ibz_init(&lattice_content);
    ibz_init(&expected_content);
    ibz_init(&bound);
    ibz_init(&norm_num);
    ibz_init(&norm_den);
    ibz_init(&norm_bound);

    parse_ideal_or_die(&secret_ideal, &cursor);
    parse_ideal_or_die(&challenge_ideal, &cursor);
    parse_ideal_or_die(&commit_ideal, &cursor);
    parse_response_or_die(&resp, &lattice_content, rtext);

    quat_lideal_inter(&challenge_secret, &challenge_ideal, &secret_ideal, &QUATALG_PINFTY);
    quat_lattice_conjugate_without_hnf(&commit_conj, &commit_ideal.lattice);
    quat_lattice_intersect(&hom_lattice, &challenge_secret.lattice, &commit_conj);
    ibz_mul(&expected_content, &challenge_secret.norm, &commit_ideal.norm);

    int ok_content = (ibz_cmp(&expected_content, &lattice_content) == 0);
    int ok_contains = quat_lattice_contains(NULL, &hom_lattice, &resp);

    ibz_pow(&bound, &ibz_const_two, SQIsign_response_length);
    ibz_sub(&bound, &bound, &ibz_const_one);
    ibz_mul(&bound, &bound, &lattice_content);
    quat_alg_norm(&norm_num, &norm_den, &resp, &QUATALG_PINFTY);
    ibz_mul(&norm_bound, &norm_den, &bound);
    int ok_bound = (ibz_cmp(&norm_num, &norm_bound) <= 0);

    printf("content=%s\n", ok_content ? "OK" : "FAIL");
    printf("contains=%s\n", ok_contains ? "OK" : "FAIL");
    printf("bound=%s\n", ok_bound ? "OK" : "FAIL");
    if (ok_content && ok_contains && ok_bound) {
        printf("OK\n");
    } else {
        printf("FAIL\n");
    }

    ibz_finalize(&norm_bound);
    ibz_finalize(&norm_den);
    ibz_finalize(&norm_num);
    ibz_finalize(&bound);
    ibz_finalize(&expected_content);
    ibz_finalize(&lattice_content);
    quat_alg_elem_finalize(&resp);
    quat_lattice_finalize(&hom_lattice);
    quat_lattice_finalize(&commit_conj);
    quat_left_ideal_finalize(&challenge_secret);
    quat_left_ideal_finalize(&commit_ideal);
    quat_left_ideal_finalize(&challenge_ideal);
    quat_left_ideal_finalize(&secret_ideal);
    free(rtext);
    free(vtext);

    return (ok_content && ok_contains && ok_bound) ? 0 : 2;
}

static int cmd_compute_vector(const char *vector_path)
{
    FILE *vf = fopen(vector_path, "rb");
    if (!vf) {
        die("open vector file failed");
    }
    char *vtext = read_all(vf);
    fclose(vf);

    char *cursor = vtext;
    quat_left_ideal_t secret_ideal;
    quat_left_ideal_t challenge_ideal;
    quat_left_ideal_t commit_ideal;
    quat_alg_elem_t resp;
    ibz_t lattice_content;

    quat_left_ideal_init(&secret_ideal);
    quat_left_ideal_init(&challenge_ideal);
    quat_left_ideal_init(&commit_ideal);
    quat_alg_elem_init(&resp);
    ibz_init(&lattice_content);

    seed_rng();
    parse_ideal_or_die(&secret_ideal, &cursor);
    parse_ideal_or_die(&challenge_ideal, &cursor);
    parse_ideal_or_die(&commit_ideal, &cursor);
    compute_response_quat_element_from_ideals(&resp,
                                             &lattice_content,
                                             &secret_ideal,
                                             &challenge_ideal,
                                             &commit_ideal);

    puts("1");
    print_ibz_line(&lattice_content);
    print_ibz_line(&resp.denom);
    for (int i = 0; i < 4; i++) {
        print_ibz_line(&resp.coord[i]);
    }

    ibz_finalize(&lattice_content);
    quat_alg_elem_finalize(&resp);
    quat_left_ideal_finalize(&commit_ideal);
    quat_left_ideal_finalize(&challenge_ideal);
    quat_left_ideal_finalize(&secret_ideal);
    free(vtext);
    return 0;
}

static int cmd_bench_vector(const char *vector_path, unsigned long count)
{
    FILE *vf = fopen(vector_path, "rb");
    if (!vf) {
        die("open vector file failed");
    }
    char *vtext = read_all(vf);
    fclose(vf);

    char *cursor = vtext;
    quat_left_ideal_t secret_ideal;
    quat_left_ideal_t challenge_ideal;
    quat_left_ideal_t commit_ideal;

    quat_left_ideal_init(&secret_ideal);
    quat_left_ideal_init(&challenge_ideal);
    quat_left_ideal_init(&commit_ideal);
    parse_ideal_or_die(&secret_ideal, &cursor);
    parse_ideal_or_die(&challenge_ideal, &cursor);
    parse_ideal_or_die(&commit_ideal, &cursor);

    seed_rng();
    cpucycles_init();
    uint64_t sum = 0;
    uint64_t min = UINT64_MAX;
    uint64_t max = 0;

    for (unsigned long i = 0; i < count; i++) {
        quat_alg_elem_t resp;
        ibz_t lattice_content;
        quat_alg_elem_init(&resp);
        ibz_init(&lattice_content);

        uint64_t t0 = cpucycles();
        compute_response_quat_element_from_ideals(&resp,
                                                 &lattice_content,
                                                 &secret_ideal,
                                                 &challenge_ideal,
                                                 &commit_ideal);
        uint64_t t1 = cpucycles();
        uint64_t dt = t1 - t0;
        sum += dt;
        if (dt < min) {
            min = dt;
        }
        if (dt > max) {
            max = dt;
        }

        ibz_finalize(&lattice_content);
        quat_alg_elem_finalize(&resp);
    }

    printf("response_core_count=%lu\n", count);
    printf("response_core_avg_cycles=%llu\n", (unsigned long long)(sum / count));
    printf("response_core_min_cycles=%llu\n", (unsigned long long)min);
    printf("response_core_max_cycles=%llu\n", (unsigned long long)max);

    quat_left_ideal_finalize(&commit_ideal);
    quat_left_ideal_finalize(&challenge_ideal);
    quat_left_ideal_finalize(&secret_ideal);
    free(vtext);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "gen") == 0) {
        return cmd_gen();
    }
    if (argc == 4 && strcmp(argv[1], "verify") == 0) {
        return cmd_verify(argv[2], argv[3]);
    }
    if (argc == 3 && strcmp(argv[1], "compute-vector") == 0) {
        return cmd_compute_vector(argv[2]);
    }
    if (argc == 4 && strcmp(argv[1], "bench-vector") == 0) {
        return cmd_bench_vector(argv[2], strtoul(argv[3], NULL, 10));
    }
    if (argc == 3 && strcmp(argv[1], "sample") == 0) {
        return cmd_sample(strtoul(argv[2], NULL, 10));
    }
    if (argc == 3 && strcmp(argv[1], "dump-responses") == 0) {
        return cmd_dump_responses(strtoul(argv[2], NULL, 10));
    }
    fprintf(stderr, "usage: %s gen > vector.txt\n", argv[0]);
    fprintf(stderr, "       %s verify vector.txt board_response.txt\n", argv[0]);
    fprintf(stderr, "       %s compute-vector vector.txt > response.txt\n", argv[0]);
    fprintf(stderr, "       %s bench-vector vector.txt N\n", argv[0]);
    fprintf(stderr, "       %s sample N > response_stats.txt\n", argv[0]);
    fprintf(stderr, "       %s dump-responses N > responses.txt\n", argv[0]);
    return 1;
}
