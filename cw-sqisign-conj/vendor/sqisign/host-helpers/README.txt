Upstream test-only helpers, used by the HOST vector generator only -- never built
into the firmware.

  helpers.c                  quat_represent_integer_even()
                             (the production quat_represent_integer() requires an
                             odd target norm; the challenge ideal needs an even
                             one, 2^CHALLENGE_BITS * m)
  random_input_generation.c  ibz_rand_interval_bits()
  quaternion_tests.h         declarations for both

Source: the-sqisign.git f417ebd, src/quaternion/ref/lvlx/test/
