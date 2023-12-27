#include "pbc/pbc.h"
#include "pbc/pbc_test.h"
#include <openssl/md5.h>
#include <openssl/sha.h>

// Number of cycle calculations
#define N_ATTR 1000

typedef struct Times_struct
{
    double mul_Z;
    double inv_Z;
    double add_Z;
    double hash;
    double add_G;
    double mul_G;
} TimesStore;

TimesStore times_Res;

int main(int argc, char **argv)
{
    pairing_t pairing;
    char param[1024];
    FILE *fp = fopen("a.param", "r");
    size_t count_param = fread(param, 1, 1024, fp);
    fclose(fp);
    pairing_init_set_buf(pairing, param, count_param);

    if (!pairing_is_symmetric(pairing))
    {
        fprintf(stderr, "only works with symmetric pairing \n");
        exit(1);
    }
    // Define elements of an algebraic structure
    element_t a, b, x, y;
    element_t mul_res_Z, inv_res_Z, add_res_Z, add_res_G, mul_res_G;

    element_init_Zr(mul_res_Z, pairing);
    element_init_Zr(inv_res_Z, pairing);
    element_init_Zr(add_res_Z, pairing);
    element_init_Zr(a, pairing);
    element_init_Zr(b, pairing);

    element_init_G1(add_res_G, pairing);
    element_init_G1(mul_res_G, pairing);
    element_init_G1(x, pairing);
    element_init_G1(y, pairing);

    element_random(a);
    element_random(b);
    element_random(x);
    element_random(y);

    double t0, t1;
    double totalmul_Z = 0.0;
    double totalinv_Z = 0.0;
    double totaladd_Z = 0.0;
    double totalhash = 0.0;
    double totaladd_G = 0.0;
    double totalmul_G = 0.0;

    for (int i = 0; i < N_ATTR; i++)
    {
        t0 = pbc_get_time();
        // mul_res_Z = a * b
        element_mul(mul_res_Z, a, b);
        t1 = pbc_get_time();
        totalmul_Z += t1 - t0;

        // evaluating time for Modular inversion operation in Z_p
        t0 = pbc_get_time();
        // Set 'inv_res_Z' to the inverse of 'a'.
        element_invert(inv_res_Z, a);
        t1 = pbc_get_time();
        totalinv_Z += t1 - t0;

        // evaluating time for Modular addition operation in Z_p
        t0 = pbc_get_time();
        // add_res_Z = a + b
        element_add(add_res_Z, a, b);
        t1 = pbc_get_time();
        totaladd_Z += t1 - t0;

        t0 = pbc_get_time();
        unsigned char *data = "123";
        unsigned char sha256_hash[SHA256_DIGEST_LENGTH];

        SHA256_CTX sha256_ctx;
        SHA256_Init(&sha256_ctx);
        SHA256_Update(&sha256_ctx, data, strlen(data));
        SHA256_Final(sha256_hash, &sha256_ctx);
        t1 = pbc_get_time();
        totalhash += t1 - t0;

        // evaluating time for Point addition in G
        t0 = pbc_get_time();
        // add_res_G = x + y
        element_add(add_res_G, x, y);
        t1 = pbc_get_time();
        totaladd_G += t1 - t0;

        // evaluating time for Point multiplication in G
        t0 = pbc_get_time();
        // mul_res_G = x * a
        element_mul_zn(mul_res_G, x, a);
        t1 = pbc_get_time();
        totalmul_G += t1 - t0;
    }

    // Divide by total number of times to calculate single time
    times_Res.mul_Z = totalmul_Z / N_ATTR;
    times_Res.inv_Z = totalinv_Z / N_ATTR;
    times_Res.add_Z = totaladd_Z / N_ATTR;
    times_Res.hash = totalhash / N_ATTR;
    times_Res.add_G = totaladd_G / N_ATTR;
    times_Res.mul_G = totalmul_G / N_ATTR;

    double tmm = times_Res.mul_Z * 1000;
    double tinv = times_Res.inv_Z * 1000;
    double tma = times_Res.add_Z * 1000;
    double th = times_Res.hash * 1000;
    double tpa = times_Res.add_G * 1000;
    double tpm = times_Res.mul_G * 1000;


    // Calculate the execution time of each cryptographic operation in Table II

    printf("                TABLE II\n");
    printf("AVERAGE RUNTIME(MS)OF CRYPTOGRAPHIC OPERATIONS\n");
    printf("+----------+-----------------------------------------------------------------------+\n");
    printf("| Notation |    Tpm    |    Tpa    |    Tmm    |    Tma    |    Tinv   |     Th    |\n");
    printf("+----------+-----------------------------------------------------------------------+\n|   Time   ");
    printf("| %.4f ms ", tpm);
    printf("| %.4f ms ", tpa);
    printf("| %.4f ms ", tmm);
    printf("| %.4f ms ", tma);
    printf("| %.4f ms ", tinv);
    printf("| %.4f ms |\n", th);
    printf("+----------+-----------------------------------------------------------------------+\n");

    // Calculate the computational overhead for invidual signatures in Table IV

    printf("                TABLE IV\n");
    printf("COMPARISON OF COMPUTATIONALCOSTS FOR INDIVIDUAL SIGNATURES\n");
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("| Scheme | Single Signature Generation(ms) | Single Signature Verification(ms) |\n");
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   19   |            %.4f ms            |            %.4f ms              |\n", tpm + 2 * tma + 2 * tmm + 2 * th, 4 * tpm + 3 * tpa + 3 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   25   |            %.4f ms            |            %.4f ms              |\n", tpm + 2 * tma + 2 * tmm + 2 * th, 3 * tpm + 2 * tpa + 2 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   27   |            %.4f ms            |            %.4f ms              |\n", tpm + tma + tmm + th, 3 * tpm + 2 * tpa + 2 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   29   |            %.4f ms            |            %.4f ms              |\n", tpm + 2 * tma + tmm + th, 3 * tpm + 3 * tpa + 2 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   30   |            %.4f ms            |            %.4f ms              |\n", tpm + tma + tmm + th, 3 * tpm + 2 * tpa + th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   32   |            %.4f ms            |            %.4f ms              |\n", tpm + 2 * tma + 2 * tmm + 2 * th, 4 * tpm + 3 * tpa + 3 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|   33   |            %.4f ms            |            %.4f ms              |\n", tpm + 2 * tma + 2 * tmm + 2 * th, 4 * tpm + 3 * tpa + 3 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");
    printf("|  Ours  |            %.4f ms            |            %.4f ms              |\n", tpm + 3 * tma + 2 * tmm + 2 * th + tinv, 4 * tpm + 3 * tpa + 3 * th);
    printf("+--------+---------------------------------+-----------------------------------+\n");

    // Calculate the computational overhead for aggregate signatures in Table V
    int n = 10;
    printf("                TABLE V\n");
    printf("COMPARISON OF COMPUTATIONAL COSTS FOR AGGREGATE SIGNATURES\n");
    printf("+--------+--------------------------------------+\n");
    printf("| Scheme | Aggregate Signature Verification(ms) |\n");
    printf("+--------+--------------------------------------+\n");
    printf("|   19   |              %.4f ms              |\n", (3 * n + 1) * tpm + (3 * n + 1) * tpa + (3 * n - 3) * tmm + 3 * n * th + (3 * n - 6) * tma);
    printf("+--------+--------------------------------------+\n");
    printf("|   25   |              %.4f ms              |\n", (2 * n + 1) * tpm + (3 * n - 1) * tpa + 2 * n * th);
    printf("+--------+--------------------------------------+\n");
    printf("|   27   |              %.4f ms              |\n", (2 * n + 2) * tpm + 2 * n * tpa + 2 * n * tmm + (n - 1) * tma + 3 * n * th);
    printf("+--------+--------------------------------------+\n");
    printf("|   29   |              %.4f ms              |\n", (n + 2) * tpm + (n + 1) * tpa + 2 * n * tmm + n * tma + 2 * n * th);
    printf("+--------+--------------------------------------+\n");
    printf("|   30   |              %.4f ms              |\n", (2 * n + 1) * tpm + (2 * n - 1) * tpa + 2 * n * tmm + 2 * n * th + (n - 1) * tma);
    printf("+--------+--------------------------------------+\n");
    printf("|   32   |              %.4f ms              |\n", (2 * n + 1) * tpm + (3 * n - 1) * tpa + 3 * n * th);
    printf("+--------+--------------------------------------+\n");
    printf("|   33   |              %.4f ms              |\n", (2 * n + 1) * tpm + 3 * n * tpa + 2 * n * th);
    printf("+--------+--------------------------------------+\n");
    printf("|  Ours  |              %.4f ms              |\n", 2 * n * tpm + (3 * n - 1) * tpa + (2 * n + 1) * th);
    printf("+--------+--------------------------------------+\n");

    // Clear variables and release occupied resources
    element_clear(a);
    element_clear(b);
    element_clear(x);
    element_clear(y);
    element_clear(mul_res_Z);
    element_clear(inv_res_Z);
    element_clear(add_res_Z);
    element_clear(add_res_G);
    element_clear(mul_res_G);
    pairing_clear(pairing);

    return 0;
}
