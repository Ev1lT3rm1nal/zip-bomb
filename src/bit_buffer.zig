const std = @import("std");

current: usize = 0,
bit_pos: usize = 0,
done: std.ArrayList(u8),

const Self = @This();

pub fn init(allocator: std.mem.Allocator) Self {
    return Self{ .done = std.ArrayList(u8).init(allocator) };
}

pub fn deinit(self: *Self) void {
    self.done.deinit();
}

pub fn push(self: *Self, x_1: u8, n_1: u8) !void {
    var x = x_1;
    var n = n_1;
    while (n >= (8 - self.bit_pos)) {
        self.current |= (x << @intCast(self.bit_pos)) & 0xff;
        x >>= @intCast(8 - self.bit_pos);
        n -= @intCast(8 - self.bit_pos);
        try self.done.append(@intCast(self.current));
        self.current = 0;
        self.bit_pos = 0;
    }
    self.current |= (x << @intCast(self.bit_pos)) & 0xff;
    self.bit_pos += n;
}

pub fn push_rev(self: *Self, x: u8, n: u8) !void {
    var mask: u8 = (@as(u8, 1) << @as(u3, @intCast(n))) >> 1;
    while (mask > 0) {
        try self.push(@intFromBool((x & mask) != 0) & 1 | 0, 1);
        mask >>= 1;
    }
}

pub fn bytes(self: *Self) ![]u8 {
    if (self.bit_pos != 0) {
        try self.done.append(@intCast(self.current));
    }
    return self.done.toOwnedSlice();
}

pub fn reset(self: *Self) void {
    self.current = 0;
    self.bit_pos = 0;
    self.done.clearAndFree();
}
