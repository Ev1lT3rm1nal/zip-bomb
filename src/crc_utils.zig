pub const std = @import("std");

pub const crc_data_type = @Vector(33, u64);

pub const crc_poly = 0xedb88320;

pub const cache_type = std.AutoArrayHashMap(u32, crc_data_type);

var crc_cache_repeated: cache_type = undefined;

pub const crc_m0: crc_data_type = blk: {
    var result: [33]u64 = undefined;
    result[0] = crc_poly;
    for (1..33) |i| {
        result[i] = 1 << (i - 1);
    }
    result[32] = 1 << 32;
    break :blk result;
};

pub const crc_mflip: crc_data_type = blk: {
    var result: [33]u64 = undefined;
    for (0..32) |i| {
        result[i] = 1 << i;
    }
    result[32] = 1 << 32 | 1;
    break :blk result;
};

pub const crc_m1 = blk: {
    @setEvalBranchQuota(2000);
    break :blk matrix_mul(crc_m0, crc_mflip);
};

pub const matrix_choose = &[_]crc_data_type{ crc_m0, crc_m1 };

pub const vector_identity: crc_data_type = @as(crc_data_type, @splat(@as(u64, 1))) << std.simd.iota(u64, 33);

pub fn vector_mul_vector(m: crc_data_type, v: u64) u64 {
    var vec = @select(u64, ((@as(crc_data_type, @splat(v)) >> comptime std.simd.iota(u64, 33)) %
        comptime @as(crc_data_type, @splat(@as(u64, 2)))) ==
        comptime @as(crc_data_type, @splat(@as(u64, 1))), m, comptime @as(crc_data_type, @splat(@as(u64, 0))));
    return @reduce(.Xor, vec);
}

pub fn matrix_mul(a: crc_data_type, b: crc_data_type) crc_data_type {
    var result: crc_data_type = undefined;
    for (0..33) |i| {
        result[i] = vector_mul_vector(a, b[i]);
    }
    return result;
}

pub fn precompute_crc_matrix(data: []u8) crc_data_type {
    var m = vector_identity;

    for (data) |b| {
        for (0..8) |shift| {
            m = matrix_mul(matrix_choose[(b >> @intCast(shift)) & 1], m);
        }
    }
    return m;
}

pub fn precompute_crc_matrix_repeated(data: []u8, n_1: usize) crc_data_type {
    var accum = precompute_crc_matrix(data);
    var n = n_1;
    var m = vector_identity;

    while (n > 0) : (n = n >> 1) {
        if (n & 1 == 1) {
            m = matrix_mul(m, accum);
        }
        accum = matrix_mul(accum, accum);
    }
    return m;
}

pub fn crc_matrix_apply(m: crc_data_type, value_1: ?u64) u64 {
    var value = value_1 orelse 0;

    return (vector_mul_vector(m, (value ^ 0xffffffff) | comptime (1 << 32)) & 0xffffffff) ^ 0xffffffff;
}

// def crc_combine(crc_prefix, crc_remainder, crc_matrix_zeroes):
//     # Undo the pre- and post-conditioning that crc_matrix_apply does, because
//     # crc_prefix and crc are already conditioned.
//     return crc_matrix_apply(crc_matrix_zeroes, crc_prefix ^ 0xffffffff) ^ 0xffffffff ^ crc_remainder

pub fn crc_combine(crc_prefix: u64, crc_remainder: u64, crc_matrix_zeroes: crc_data_type) u64 {
    return crc_matrix_apply(crc_matrix_zeroes, crc_prefix ^ 0xffffffff) ^ 0xffffffff ^ crc_remainder;
}

pub fn hash_data(data: []u8, n: usize) u32 {
    var num: [@bitSizeOf(usize) / 8]u8 = @bitCast(n);
    var hasher = std.hash.Crc32.init();
    hasher.update(data);
    hasher.update("@");
    hasher.update(&num);
    return hasher.final();
}

pub fn memo_crc_matrix_repeated(data: []u8, n: usize) crc_data_type {
    var key = hash_data(data, n);
    var value = crc_cache_repeated.get(key);
    var result: crc_data_type = undefined;
    if (value) |final| {
        result = final;
    } else {
        result = precompute_crc_matrix_repeated(data, n);
        _ = crc_cache_repeated.put(key, result) catch {};
    }
    return result;
}

pub fn init_cache(allocator: std.mem.Allocator) void {
    crc_cache_repeated = cache_type.init(allocator);
}

pub fn deinit_cache() void {
    crc_cache_repeated.deinit();
}
