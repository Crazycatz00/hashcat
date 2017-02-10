/**
 * Author......: Crazycatz00
 * License.....: MIT
 */

#define _CRC32_
//#define NEW_SIMD_CODE

#include "inc_vendor.cl"
#include "inc_types.cl"
#include "inc_common.cl"
#include "inc_rp.h"
#include "inc_rp.cl"
#include "inc_simd.cl"

#if VECT_SIZE == 1
#define CONVERTX(type)
#else
#define CONVERTX3(type, width) convert_ ## type ## width
#define CONVERTX2(type, width) CONVERTX3(type, width)
#define CONVERTX(type) CONVERTX2(type, VECT_SIZE)
#endif

static u32x round_kh2hp(u32x a, u32x c)
{
  a ^= c << 24;
  for (u32 i = 0; i < 8; i += 1)
  {
    a = (a & 0x80000000) != 0 ? (a << 1) ^ 0x04c11db7 : a << 1;
  }

  return a;
}

static u32x kh2hp(const u32x w[16], const u32 pw_len)
{
  u32x a = ~0;

  if (pw_len >= 1) a = round_kh2hp(a, w[0] >>  0);
  if (pw_len >= 2) a = round_kh2hp(a, w[0] >>  8);
  if (pw_len >= 3) a = round_kh2hp(a, w[0] >> 16);
  if (pw_len >= 4) a = round_kh2hp(a, w[0] >> 24);

  if (pw_len >= 5) a = round_kh2hp(a, w[1] >>  0);
  if (pw_len >= 6) a = round_kh2hp(a, w[1] >>  8);
  if (pw_len >= 7) a = round_kh2hp(a, w[1] >> 16);
  if (pw_len >= 8) a = round_kh2hp(a, w[1] >> 24);

  for (u32 i = 8, j = 2; i < pw_len; i += 4, j += 1)
  {
    if (pw_len >= (i + 1)) a = round_kh2hp(a, w[j] >>  0);
    if (pw_len >= (i + 2)) a = round_kh2hp(a, w[j] >>  8);
    if (pw_len >= (i + 3)) a = round_kh2hp(a, w[j] >> 16);
    if (pw_len >= (i + 4)) a = round_kh2hp(a, w[j] >> 24);
  }

  return ~a;
}

static u32x round_kh2hs(u32x a, u32x c)
{
  a ^= c << 8;
  for (u32 i = 0; i < 8; i += 1)
  {
    a = (a & 0x8000) != 0 ? (a << 1) ^ 0x1021 : a << 1;
  }

  return a;
}

static u16x kh2hs(const u32x w[16], const u32 pw_len)
{
  u32x a = ~0;

  if (pw_len >= 1)
  {
    for (u32 j = (pw_len - 1) / 4, i = j * 4; ; i -= 4, j -= 1)
    {
      if (pw_len >= (i + 4)) a = round_kh2hs(a, w[j] >> 24);
      if (pw_len >= (i + 3)) a = round_kh2hs(a, w[j] >> 16);
      if (pw_len >= (i + 2)) a = round_kh2hs(a, w[j] >>  8);
      if (pw_len >= (i + 1)) a = round_kh2hs(a, w[j] >>  0);
      if (i == 0) break;
    }
  }

  return CONVERTX(ushort)(~a);
}

__kernel void m01313_m04(__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  const u32 gid = get_global_id(0);
  if (gid >= gid_max) return;
  const u32 lid = get_local_id(0);  

  u32 pw_buf0[4] =
  {
    pws[gid].i[0],
    pws[gid].i[1],
    pws[gid].i[2],
    pws[gid].i[3]
  };
  u32 pw_buf1[4] =
  {
    pws[gid].i[4],
    pws[gid].i[5],
    pws[gid].i[6],
    pws[gid].i[7]
  };

  const u32 pw_len = pws[gid].pw_len;
  const u32x z = 0;

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };

    const u32x out_len = apply_rules_vect(pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    const u32x w_t[16] =
    {
	  w0[0],
	  w0[1],
	  w0[2],
	  w0[3],
	  w1[0],
	  w1[1],
	  w1[2],
	  w1[3],
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
    };
    const u32x hashP = kh2hp(w_t, pw_len);
    const u16x hashS = kh2hs(w_t, pw_len);

    COMPARE_M_SIMD(hashP, hashS, z, z);
  }
}

__kernel void m01313_m08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m01313_m16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m01313_s04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
  const u32 gid = get_global_id(0);
  if (gid >= gid_max) return;
  const u32 lid = get_local_id(0);
  
  u32 pw_buf0[4] =
  {
    pws[gid].i[0],
    pws[gid].i[1],
    pws[gid].i[2],
    pws[gid].i[3]
  };
  u32 pw_buf1[4] =
  {
    pws[gid].i[4],
    pws[gid].i[5],
    pws[gid].i[6],
    pws[gid].i[7]
  };

  const u32 search[4] =
  {
    digests_buf[digests_offset].digest_buf[DGST_R0],
    digests_buf[digests_offset].digest_buf[DGST_R1],
    0,
    0
  };
  const u32 pw_len = pws[gid].pw_len;
  const u32x z = 0;

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    u32x w0[4] = { 0 };
    u32x w1[4] = { 0 };

    const u32x out_len = apply_rules_vect(pw_buf0, pw_buf1, pw_len, rules_buf, il_pos, w0, w1);

    const u32x w_t[16] =
    {
	  w0[0],
	  w0[1],
	  w0[2],
	  w0[3],
	  w1[0],
	  w1[1],
	  w1[2],
	  w1[3],
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
	  0,
    };
    const u32x hashP = kh2hp(w_t, pw_len);
    if (MATCHES_NONE_VS(hashP, search[0])) continue;
    const u16x hashS = kh2hs(w_t, pw_len);

    COMPARE_S_SIMD(hashP, hashS, z, z);
  }
}

__kernel void m01313_s08 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}

__kernel void m01313_s16 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
{
}
