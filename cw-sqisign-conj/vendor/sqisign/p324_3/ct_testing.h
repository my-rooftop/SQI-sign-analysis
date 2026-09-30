// SPDX-License-Identifier: Apache-2.0

#ifndef CT_TESTING_H
#define CT_TESTING_H

// Constant-time (ctgrind/valgrind) taint-tracking annotations.
//
// CT_TESTING_MAKE_PUBLIC(addr, len)  declassify: tell valgrind the region is defined/public (e.g. a secret-derived
//                                                value argued to be safe to branch on).
// CT_TESTING_MAKE_SECRET(addr, len)  poison:     tell valgrind the region is undefined/secret (used at the RNG roots
//                                                and to taint serialized secret keys).
//
// The macros are active only when ENABLE_CT_TESTING is defined (a global compile definition set by the  build). When it
// is not defined (the default), both expand to a no-op that still evaluates (casts to void) each argument, so a
// variable whose only use is one of these annotations does not trip -Wunused-variable under -Werror.

#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>

#define CT_TESTING_MAKE_PUBLIC(addr, len) VALGRIND_MAKE_MEM_DEFINED((addr), (len))
#define CT_TESTING_MAKE_SECRET(addr, len) VALGRIND_MAKE_MEM_UNDEFINED((addr), (len))
#else
#define CT_TESTING_MAKE_PUBLIC(addr, len) ((void)(addr), (void)(len))
#define CT_TESTING_MAKE_SECRET(addr, len) ((void)(addr), (void)(len))
#endif

#endif /* CT_TESTING_H */
