#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

gmp_randstate_t randomGeneratorState;
unsigned long randomSeed;

struct PubKey
{
    mpz_t n;
    mpz_t nSquared;
    mpz_t g;
};

struct PrivKey
{
    mpz_t p;
    mpz_t q;
    mpz_t pMinusOne;
    mpz_t qMinusOne;
    mpz_t pSquared;
    mpz_t qSquared;
    mpz_t hp;
    mpz_t hq;
    mpz_t qTimesQInvModP;
    mpz_t pTimesPInvModQ;
    mpz_t posNegBoundary;
};

void initPubKey(struct PubKey *pubKey)
{
    mpz_init(pubKey->n);
    mpz_init(pubKey->nSquared);
    mpz_init(pubKey->g);
}

void freePubKey(struct PubKey *pubKey)
{
    mpz_clear(pubKey->n);
    mpz_clear(pubKey->nSquared);
    mpz_clear(pubKey->g);
}

void initPrivKey(struct PrivKey *privKey)
{
    mpz_init(privKey->p);
    mpz_init(privKey->q);
    mpz_init(privKey->pMinusOne);
    mpz_init(privKey->qMinusOne);
    mpz_init(privKey->pSquared);
    mpz_init(privKey->qSquared);
    mpz_init(privKey->hp);
    mpz_init(privKey->hq);
    mpz_init(privKey->qTimesQInvModP);
    mpz_init(privKey->pTimesPInvModQ);
    mpz_init(privKey->posNegBoundary);
}

void freePrivKey(struct PrivKey *privKey)
{
    mpz_clear(privKey->p);
    mpz_clear(privKey->q);
    mpz_clear(privKey->pMinusOne);
    mpz_clear(privKey->qMinusOne);
    mpz_clear(privKey->pSquared);
    mpz_clear(privKey->qSquared);
    mpz_clear(privKey->hp);
    mpz_clear(privKey->hq);
    mpz_clear(privKey->qTimesQInvModP);
    mpz_clear(privKey->pTimesPInvModQ);
    mpz_clear(privKey->posNegBoundary);
}

unsigned long getRandomSeed()
{
    FILE *urandom;
    unsigned long randomNumber;

    urandom = fopen("/dev/urandom", "r");

    if (urandom == NULL)
    {
        fprintf(stderr, "Cannot open /dev/urandom!\n");
        exit(1);
    }

    fread(&randomNumber, sizeof(randomNumber), 1, urandom);

    return randomNumber;
}

void getRandomPrime(mpz_t output, unsigned long primeLen)
{
    do
    {
        //generate a random number in the interval [0, 2^(numberOfBits - 1))
        mpz_urandomb(output, randomGeneratorState, primeLen - 1);

        //shift number to the interval [2^(numberOfBits - 1), 2^numberOfBits)
        mpz_setbit(output, primeLen - 1);
    }
    while (!mpz_probab_prime_p(output, 10));
}

void print(mpz_t input)
{
    unsigned int base = 10;

    //get a pointer to GMP's internal memory deallocator function
    void (*deallocator)(void *, size_t);
    mp_get_memory_functions(NULL, NULL, &deallocator);

    //get the string representation of input
    char *data = mpz_get_str(NULL, base, input);

    printf("%s\n", data);

    (*deallocator)((void *)data, strlen(data));
}

//Computes @f$ L(u) = \frac{u - 1}{d} @f$
void L(mpz_t output, mpz_t input, mpz_t d)
{
    mpz_sub_ui(output, input, 1);

    mpz_tdiv_q(output, output, d);
}

void generateKeys(struct PubKey *pubKey, struct PrivKey *privKey)
{
    unsigned long keyLen = 64;
    unsigned long primeLen = keyLen / 2;

    mpz_t tmp;
    mpz_init(tmp);

    do
    {
        getRandomPrime(privKey->p, primeLen);
        getRandomPrime(privKey->q, primeLen);

        while (!mpz_cmp(privKey->p, privKey->q))
        {
            getRandomPrime(privKey->p, primeLen);
        }

        mpz_mul(pubKey->n, privKey->p, privKey->q);
    }
    while (mpz_sizeinbase(pubKey->n, 2) != keyLen);

    mpz_mul(pubKey->nSquared, pubKey->n, pubKey->n);

    mpz_add_ui(pubKey->g, pubKey->n, 1);

    mpz_sub_ui(privKey->pMinusOne, privKey->p, 1);
    mpz_sub_ui(privKey->qMinusOne, privKey->q, 1);

    mpz_mul(privKey->pSquared, privKey->p, privKey->p);
    mpz_mul(privKey->qSquared, privKey->q, privKey->q);

    //compute hp
    mpz_powm(tmp, pubKey->g, privKey->pMinusOne, privKey->pSquared);
    L(privKey->hp, tmp, privKey->p);
    mpz_invert(privKey->hp, privKey->hp, privKey->p);

    //compute hq
    mpz_powm(tmp, pubKey->g, privKey->qMinusOne, privKey->qSquared);
    L(privKey->hq, tmp, privKey->q);
    mpz_invert(privKey->hq, privKey->hq, privKey->q);

    //precomputations
    mpz_invert(tmp, privKey->p, privKey->q);
    mpz_mul(privKey->pTimesPInvModQ, privKey->p, tmp);
    mpz_invert(tmp, privKey->q, privKey->p);
    mpz_mul(privKey->qTimesQInvModP, privKey->q, tmp);

    mpz_tdiv_q_ui(privKey->posNegBoundary, pubKey->n, 2);

    mpz_clear(tmp);

    printf("n: ");
    print(pubKey->n);
    printf("Positive / negative boundary: ");
    print(privKey->posNegBoundary);
}

void encrypt_ul(mpz_t output, long input, struct PubKey *pubKey)
{
    mpz_t tmp;
    mpz_init(tmp);

    if (input < 0)
    {
        mpz_sub_ui(tmp, pubKey->n, -input);
    }
    else
    {
        mpz_set_ui(tmp, input);
    }

    mpz_mul(output, pubKey->n, tmp);
    mpz_add_ui(output, output, 1);

    //mpz_mod(output, output, pubKey->nSquared);//the first step of decryption is % n^2 anyway...

    mpz_clear(tmp);
}

void decrypt(mpz_t output, mpz_t input, struct PubKey *pubKey, struct PrivKey *privKey)
{
    mpz_t mp, mq, tmp, tmp2;
    mpz_init(mp);
    mpz_init(mq);
    mpz_init(tmp);
    mpz_init(tmp2);

    mpz_powm(tmp, input, privKey->pMinusOne, privKey->pSquared);
    L(tmp, tmp, privKey->p);
    mpz_mul(tmp, tmp, privKey->hp);
    mpz_mod(mp, tmp, privKey->p);

    mpz_powm(tmp, input, privKey->qMinusOne, privKey->qSquared);
    L(tmp, tmp, privKey->q);
    mpz_mul(tmp, tmp, privKey->hq);
    mpz_mod(mq, tmp, privKey->q);

    mpz_mul(tmp, mp, privKey->qTimesQInvModP);
    mpz_mul(tmp2, mq, privKey->pTimesPInvModQ);
    mpz_add(tmp, tmp, tmp2);
    mpz_mod(output, tmp, pubKey->n);

    if (mpz_cmp(output, privKey->posNegBoundary) > 0)
    {
        mpz_sub(output, output, pubKey->n);
    }

    mpz_clear(mp);
    mpz_clear(mq);
    mpz_clear(tmp);
    mpz_clear(tmp2);
}

void testHomomorphicSubtraction(long lhs, long rhs, struct PubKey *pubKey, struct PrivKey *privKey)
{
    //encryptions
    mpz_t encLhs, encRhs;
    mpz_init(encLhs);
    mpz_init(encRhs);

    //descryptions
    mpz_t decLhs, decRhs;
    mpz_init(decLhs);
    mpz_init(decRhs);

    //homomorphic operation result
    mpz_t encSub, decSub;
    mpz_init(encSub);
    mpz_init(decSub);

    mpz_t tmp;
    mpz_init(tmp);

    //encrypt
    encrypt_ul(encLhs, lhs, pubKey);
    encrypt_ul(encRhs, rhs, pubKey);

    //decrypt
    decrypt(decLhs, encLhs, pubKey, privKey);
    decrypt(decRhs, encRhs, pubKey, privKey);

    printf("Decrypted lhs: ");
    print(decLhs);
    printf("Decrypted rhs: ");
    print(decRhs);

    //homomorphic subtraction
    mpz_invert(tmp, encRhs, pubKey->nSquared);
    mpz_mul(encSub, encLhs, tmp);
    //mpz_mod(sub, sub, nSquared);//the first step of decryption is % n^2 anyway...

    decrypt(decSub, encSub, pubKey, privKey);

    printf("Decrypted subtraction: ");
    print(decSub);

    mpz_clear(encLhs);
    mpz_clear(encRhs);
    mpz_clear(decLhs);
    mpz_clear(decRhs);
    mpz_clear(encSub);
    mpz_clear(decSub);

    mpz_clear(tmp);
}

int main(int argc, char *argv[])
{
    randomSeed = getRandomSeed();
    gmp_randinit_default(randomGeneratorState);
    gmp_randseed_ui(randomGeneratorState, randomSeed);

    struct PubKey pubKey;
    struct PrivKey privKey;
    initPubKey(&pubKey);
    initPrivKey(&privKey);
    generateKeys(&pubKey, &privKey);

    long lhs = strtol(argv[1], NULL, 10);
    long rhs = strtol(argv[2], NULL, 10);

    printf("Testing: %ld - %ld\n", lhs, rhs);

    testHomomorphicSubtraction(lhs, rhs, &pubKey, &privKey);

    freePubKey(&pubKey);
    freePrivKey(&privKey);

    return 0;
}