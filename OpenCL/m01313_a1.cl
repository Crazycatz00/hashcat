/**
 * Author......: Crazycatz00
 * License.....: MIT
 */

//#define NEW_SIMD_CODE

#include "inc_vendor.cl"
#include "inc_types.cl"
#include "inc_common.cl"
#include "inc_simd.cl"

#if VECT_SIZE == 1
#define CONVERTX(type)
#else
#define CONVERTX3(type, width) convert_ ## type ## width
#define CONVERTX2(type, width) CONVERTX3(type, width)
#define CONVERTX(type) CONVERTX2(type, VECT_SIZE)
#endif

u32x round_kh2hp(u32x a, u32x c)
{
  a ^= c << 24;
  for (u32 i = 0; i < 8; i += 1)
  {
    a = (a & 0x80000000) != 0 ? (a << 1) ^ 0x04c11db7 : a << 1;
  }

  return a;
}

u32x kh2hp(const u32x w[16], const u32 pw_len)
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

u32x round_kh2hs(u32x a, u32x c)
{
  a ^= c << 8;
  for (u32 i = 0; i < 8; i += 1)
  {
    a = (a & 0x8000) != 0 ? (a << 1) ^ 0x1021 : a << 1;
  }

  return a;
}

u16x kh2hs(const u32x w[16], const u32 pw_len)
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

__kernel void m01313_m04 (__global pw_t *pws, __global const kernel_rule_t *rules_buf, __global const comb_t *combs_buf, __global const bf_t *bfs_buf, __global void *tmps, __global void *hooks, __global const u32 *bitmaps_buf_s1_a, __global const u32 *bitmaps_buf_s1_b, __global const u32 *bitmaps_buf_s1_c, __global const u32 *bitmaps_buf_s1_d, __global const u32 *bitmaps_buf_s2_a, __global const u32 *bitmaps_buf_s2_b, __global const u32 *bitmaps_buf_s2_c, __global const u32 *bitmaps_buf_s2_d, __global plain_t *plains_buf, __global const digest_t *digests_buf, __global u32 *hashes_shown, __global const salt_t *salt_bufs, __global const void *esalt_bufs, __global u32 *d_return_buf, __global u32 *d_scryptV0_buf, __global u32 *d_scryptV1_buf, __global u32 *d_scryptV2_buf, __global u32 *d_scryptV3_buf, const u32 bitmap_mask, const u32 bitmap_shift1, const u32 bitmap_shift2, const u32 salt_pos, const u32 loop_pos, const u32 loop_cnt, const u32 il_cnt, const u32 digests_cnt, const u32 digests_offset, const u32 combs_mode, const u32 gid_max)
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

  const u32 pw_l_len = pws[gid].pw_len;
  const u32x z = 0;

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    const u32x pw_r_len = pwlenx_create_combt(combs_buf, il_pos);
    const u32x pw_len = pw_l_len + pw_r_len;

    u32x wordl0[4] =
	{
		pw_buf0[0],
		pw_buf0[1],
		pw_buf0[2],
		pw_buf0[3]
	};
    u32x wordl1[4] =
	{
		pw_buf1[0],
		pw_buf1[1],
		pw_buf1[2],
		pw_buf1[3]
	};
    u32x wordl2[4] = { 0 };
    u32x wordl3[4] = { 0 };

    u32x wordr0[4] =
	{
		ix_create_combt(combs_buf, il_pos, 0),
		ix_create_combt(combs_buf, il_pos, 1),
		ix_create_combt(combs_buf, il_pos, 2),
		ix_create_combt(combs_buf, il_pos, 3)
	};
    u32x wordr1[4] =
	{
		ix_create_combt(combs_buf, il_pos, 4),
		ix_create_combt(combs_buf, il_pos, 5),
		ix_create_combt(combs_buf, il_pos, 6),
		ix_create_combt(combs_buf, il_pos, 7)
	};
    u32x wordr2[4] = { 0 };
    u32x wordr3[4] = { 0 };

    if (combs_mode == COMBINATOR_MODE_BASE_LEFT)
    {
      switch_buffer_by_offset_le_VV(wordr0, wordr1, wordr2, wordr3, pw_l_len);
    }
    else
    {
      switch_buffer_by_offset_le_VV(wordl0, wordl1, wordl2, wordl3, pw_r_len);
    }

    const u32x w_t[16] =
    {
	  wordl0[0] | wordr0[0],
	  wordl0[1] | wordr0[1],
	  wordl0[2] | wordr0[2],
	  wordl0[3] | wordr0[3],
	  wordl1[0] | wordr1[0],
	  wordl1[1] | wordr1[1],
	  wordl1[2] | wordr1[2],
	  wordl1[3] | wordr1[3],
	  wordl2[0] | wordr2[0],
	  wordl2[1] | wordr2[1],
	  wordl2[2] | wordr2[2],
	  wordl2[3] | wordr2[3],
	  wordl3[0] | wordr3[0],
	  wordl3[1] | wordr3[1],
	  wordl3[2] | wordr3[2],
	  wordl3[3] | wordr3[3]
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
  const u32 pw_l_len = pws[gid].pw_len;
  const u32x z = 0;

  for (u32 il_pos = 0; il_pos < il_cnt; il_pos += VECT_SIZE)
  {
    const u32x pw_r_len = pwlenx_create_combt(combs_buf, il_pos);
    const u32x pw_len = pw_l_len + pw_r_len;

    u32x wordl0[4] =
	{
		pw_buf0[0],
		pw_buf0[1],
		pw_buf0[2],
		pw_buf0[3]
	};
    u32x wordl1[4] =
	{
		pw_buf1[0],
		pw_buf1[1],
		pw_buf1[2],
		pw_buf1[3]
	};
    u32x wordl2[4] = { 0 };
    u32x wordl3[4] = { 0 };
	
    u32x wordr0[4] =
	{
		ix_create_combt(combs_buf, il_pos, 0),
		ix_create_combt(combs_buf, il_pos, 1),
		ix_create_combt(combs_buf, il_pos, 2),
		ix_create_combt(combs_buf, il_pos, 3)
	};
    u32x wordr1[4] =
	{
		ix_create_combt(combs_buf, il_pos, 4),
		ix_create_combt(combs_buf, il_pos, 5),
		ix_create_combt(combs_buf, il_pos, 6),
		ix_create_combt(combs_buf, il_pos, 7)
	};
    u32x wordr2[4] = { 0 };
    u32x wordr3[4] = { 0 };

    if (combs_mode == COMBINATOR_MODE_BASE_LEFT)
    {
      switch_buffer_by_offset_le_VV(wordr0, wordr1, wordr2, wordr3, pw_l_len);
    }
    else
    {
      switch_buffer_by_offset_le_VV(wordl0, wordl1, wordl2, wordl3, pw_r_len);
    }

    const u32x w_t[16] =
    {
	  wordl0[0] | wordr0[0],
	  wordl0[1] | wordr0[1],
	  wordl0[2] | wordr0[2],
	  wordl0[3] | wordr0[3],
	  wordl1[0] | wordr1[0],
	  wordl1[1] | wordr1[1],
	  wordl1[2] | wordr1[2],
	  wordl1[3] | wordr1[3],
	  wordl2[0] | wordr2[0],
	  wordl2[1] | wordr2[1],
	  wordl2[2] | wordr2[2],
	  wordl2[3] | wordr2[3],
	  wordl3[0] | wordr3[0],
	  wordl3[1] | wordr3[1],
	  wordl3[2] | wordr3[2],
	  wordl3[3] | wordr3[3]
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
