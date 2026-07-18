#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SQISIGN_RAM_BUDGET
#define SQISIGN_RAM_BUDGET (100u * 1024u)
#endif

#ifndef SQISIGN_RAM_STATIC
#define SQISIGN_RAM_STATIC 13240u
#endif

static uintptr_t stack_base;
static size_t heap_current;
static size_t heap_peak;
static size_t stack_peak;
static int in_wrap;

void *__real_malloc(size_t size);
void *__real_calloc(size_t nmemb, size_t size);
void *__real_realloc(void *ptr, size_t size);
void __real_free(void *ptr);
void __wrap_free(void *ptr);

static void *gmp_alloc(size_t size);
static void *gmp_realloc(void *ptr, size_t old_size, size_t new_size);
static void gmp_free(void *ptr, size_t size);
void __gmp_set_memory_functions(void *(*alloc_func)(size_t),
                                void *(*realloc_func)(void *, size_t, size_t),
                                void (*free_func)(void *, size_t)) __attribute__((weak));
void mp_set_memory_functions(void *(*alloc_func)(size_t),
                             void *(*realloc_func)(void *, size_t, size_t),
                             void (*free_func)(void *, size_t)) __attribute__((weak));
static void ram_probe_init(void) __attribute__((constructor));

static void
ram_probe_init(void)
{
    volatile unsigned char marker;
    stack_base = (uintptr_t)&marker;
    if (__gmp_set_memory_functions) {
        __gmp_set_memory_functions(gmp_alloc, gmp_realloc, gmp_free);
    } else if (mp_set_memory_functions) {
        mp_set_memory_functions(gmp_alloc, gmp_realloc, gmp_free);
    }
}

static size_t
stack_used_now(void)
{
    volatile unsigned char marker;
    uintptr_t sp = (uintptr_t)&marker;

    if (stack_base >= sp) {
        return (size_t)(stack_base - sp);
    }
    return (size_t)(sp - stack_base);
}

static size_t
checked_mul(size_t a, size_t b)
{
    if (a != 0 && b > ((size_t)-1) / a) {
        return (size_t)-1;
    }
    return a * b;
}

static void *
track_alloc(size_t size)
{
    size_t total = size + sizeof(size_t);
    size_t *raw;

    if (total < size) {
        return NULL;
    }
    raw = __real_malloc(total);
    if (!raw) {
        return NULL;
    }
    *raw = size;
    heap_current += size;
    if (heap_current > heap_peak) {
        heap_peak = heap_current;
    }
    return raw + 1;
}

static void *
gmp_alloc(size_t size)
{
    return track_alloc(size);
}

static void *
gmp_realloc(void *ptr, size_t old_size, size_t new_size)
{
    void *ret;

    (void)old_size;
    if (!ptr) {
        return track_alloc(new_size);
    }
    if (new_size == 0) {
        __wrap_free(ptr);
        return NULL;
    }
    ret = track_alloc(new_size);
    if (ret) {
        size_t tracked_old_size = ((size_t *)ptr)[-1];
        memcpy(ret, ptr, tracked_old_size < new_size ? tracked_old_size : new_size);
        heap_current -= tracked_old_size;
        __real_free(((size_t *)ptr) - 1);
    }
    return ret;
}

static void
gmp_free(void *ptr, size_t size)
{
    (void)size;
    __wrap_free(ptr);
}

void *
__wrap_malloc(size_t size)
{
    if (in_wrap) {
        return __real_malloc(size);
    }
    in_wrap = 1;
    void *ret = track_alloc(size);
    in_wrap = 0;
    return ret;
}

void *
__wrap_calloc(size_t nmemb, size_t size)
{
    size_t total;
    void *ret;

    if (in_wrap) {
        return __real_calloc(nmemb, size);
    }
    total = checked_mul(nmemb, size);
    if (total == (size_t)-1) {
        return NULL;
    }
    in_wrap = 1;
    ret = track_alloc(total);
    if (ret) {
        memset(ret, 0, total);
    }
    in_wrap = 0;
    return ret;
}

void *
__wrap_realloc(void *ptr, size_t size)
{
    size_t old_size = 0;
    void *ret;

    if (in_wrap) {
        return __real_realloc(ptr, size);
    }
    if (!ptr) {
        return __wrap_malloc(size);
    }
    if (size == 0) {
        __wrap_free(ptr);
        return NULL;
    }

    old_size = ((size_t *)ptr)[-1];
    in_wrap = 1;
    ret = track_alloc(size);
    if (ret) {
        memcpy(ret, ptr, old_size < size ? old_size : size);
        heap_current -= old_size;
        __real_free(((size_t *)ptr) - 1);
    }
    in_wrap = 0;
    return ret;
}

void
__wrap_free(void *ptr)
{
    size_t *raw;

    if (in_wrap) {
        __real_free(ptr);
        return;
    }
    if (!ptr) {
        return;
    }
    raw = ((size_t *)ptr) - 1;
    heap_current -= *raw;
    __real_free(raw);
}

void
sqisign_cw_progress(unsigned char tag)
{
    size_t stack_used = stack_used_now();
    size_t total_current;
    size_t total_peakish;

    if (stack_used > stack_peak) {
        stack_peak = stack_used;
    }
    total_current = SQISIGN_RAM_STATIC + heap_current + stack_used;
    total_peakish = SQISIGN_RAM_STATIC + heap_peak + stack_peak;

    printf("RAM tag=%c static=%u stack=%zu stack_peak=%zu heap=%zu heap_peak=%zu total=%zu peakish=%zu budget=%u%s\n",
           tag,
           (unsigned)SQISIGN_RAM_STATIC,
           stack_used,
           stack_peak,
           heap_current,
           heap_peak,
           total_current,
           total_peakish,
           (unsigned)SQISIGN_RAM_BUDGET,
           total_peakish > SQISIGN_RAM_BUDGET ? " OVER" : "");
    fflush(stdout);
}
