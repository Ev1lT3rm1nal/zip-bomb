const std = @import("std");
const native_endian = @import("builtin").target.cpu.arch.endian();
const crc_u = @import("crc_utils.zig");

const zip_version: u16 = 20;
const zip64_version: u16 = 45;

const zip_bzip2_version: u16 = 46;

pub const compression_method_deflate: u16 = 8;
pub const compression_method_bzip: u16 = 12;

const mod_date: u16 = 0x0548;
const mod_time: u16 = 0x6ca0;

pub fn zip_version_num(method: u16, zip64: bool) u16 {
    if (zip64) {
        return zip64_version;
    } else if (method == compression_method_bzip) {
        return zip_bzip2_version;
    } else {
        return zip_version;
    }
}

pub const Quote = packed struct {
    flag: u8,
    num_quote: u16,
    deflate: u16,

    pub fn new(flag: u8, num_quote: u16, deflate: u16) Quote {
        return .{
            .flag = flag,
            .num_quote = num_quote,
            .deflate = deflate,
        };
    }

    pub fn toBytes(self: Quote, allocator: std.mem.Allocator) ![]u8 {
        var bytes: [@bitSizeOf(@This()) / 8]u8 = @bitCast(self);
        if (comptime native_endian != .little) {
            std.mem.reverse(u8, &bytes);
        }
        const result = try allocator.alloc(u8, @bitSizeOf(@This()) / 8);
        @memcpy(result, &bytes);
        return result;
    }
};

const ExtraDataMin = packed struct {
    tag: u16,
    size: u16,
    uncompressed_size: u64,
};

const ExtraData = packed struct {
    min: ExtraDataMin,
    compressed_size: u64,
};

const ExtraLength = packed struct {
    tag: u16,
    size: u16,
};

const FileHeaderData = packed struct {
    signature: u32,
    version: u16,
    flags: u16,
    compression_method: u16,
    mod_time: u16,
    mod_date: u16,
    crc: u32,
    compressed_size: u32,
    uncompressed_size: u32,
    filename_length: u16,
    extra_length: u16,
};

const CentralDirectoryData = packed struct {
    signature: u32,
    version_by: u16,
    version: u16,
    flags: u16,
    compression_method: u16,
    mod_time: u16,
    mod_date: u16,
    crc: u32,
    compressed_size: u32,
    uncompressed_size: u32,
    filename_length: u16,
    extra_length: u16,
    comment_size: u16,
    start_disk: u16,
    internal_attributes: u16,
    external_attributes: u32,
    offset: u32,
};

pub const LocalFileHeader = struct {
    compressed_size: u64,
    uncompressed_size: u64,
    crc: u32,
    filename: []u8,
    compression_method: u16,
    extra_tag: ?u16,
    extra_length_excess: u16,

    pub fn init(
        compressed_size: u64,
        uncompressed_size: u64,
        crc: u32,
        filename: []u8,
        compression_method: ?u16,
        extra_tag: ?u16,
        extra_length_excess: ?u16,
    ) LocalFileHeader {
        const self = .{
            .compressed_size = compressed_size,
            .uncompressed_size = uncompressed_size,
            .crc = crc,
            .filename = filename,
            .compression_method = compression_method orelse compression_method_deflate,
            .extra_tag = extra_tag,
            .extra_length_excess = extra_length_excess orelse 0,
        };

        return self;
    }

    pub fn deinit(self: *LocalFileHeader, allocator: std.mem.Allocator) void {
        allocator.free(self.filename);
    }

    pub fn serialize(self: *LocalFileHeader, zip64: bool, allocator: std.mem.Allocator) ![]u8 {
        var extra: ?ExtraData = null;
        const extra_data_size = @bitSizeOf(ExtraData) / 8;
        const extra_length_size = @bitSizeOf(ExtraLength) / 8;
        var extra_length: ?ExtraLength = null;
        if (zip64) {
            extra = ExtraData{
                .min = .{
                    .tag = 0x0001,
                    .size = 16,
                    .uncompressed_size = self.uncompressed_size,
                },
                .compressed_size = self.compressed_size,
            };
        }

        if (self.extra_length_excess > 0) {
            extra_length = ExtraLength{
                .tag = self.extra_tag.?,
                .size = self.extra_length_excess,
            };
        }
        const compressed_size: u32 = @max(if (zip64) @as(u32, 0xffffffff) else @as(u32, 0), @as(u32, @truncate(self.compressed_size)));
        // std.debug.print("size {d}\n", .{self.uncompressed_size});
        const uncompressed_size: u32 = @max(if (zip64) @as(u32, 0xffffffff) else @as(u32, 0), @as(u32, @truncate(self.uncompressed_size)));

        const header: FileHeaderData = .{
            .signature = 0x04034b50,
            .version = zip_version_num(self.compression_method, zip64),
            .flags = 0,
            .compression_method = self.compression_method,
            .mod_time = mod_time,
            .mod_date = mod_date,
            .crc = self.crc,
            .compressed_size = compressed_size,
            .uncompressed_size = uncompressed_size,
            .filename_length = @as(u16, @intCast(self.filename.len)),
            .extra_length = (if (extra != null) @as(u16, extra_data_size) else 0) + (if (extra_length != null) @as(u16, extra_length_size) else 0),
        };
        var header_bytes: [@bitSizeOf(FileHeaderData) / 8]u8 = @bitCast(header);
        var extra_bytes: [extra_data_size]u8 = if (extra != null) @bitCast(extra.?) else [_]u8{0} ** (extra_data_size);
        var extra_length_bytes: [extra_length_size]u8 = if (extra_length != null) @bitCast(extra_length.?) else [_]u8{0} ** (extra_length_size);
        const extra_length_bytes_size: usize = if (extra_length != null) extra_length_size else 0;
        const extra_data_bytes_size: usize = if (extra != null) extra_data_size else 0;
        if (comptime native_endian != .little) {
            std.mem.reverse(u8, &header_bytes);

            if (zip64) {
                std.mem.reverse(u8, &extra_bytes);
            }
            if (self.extra_length_excess > 0) {
                std.mem.reverse(u8, &extra_length_bytes);
            }
        }
        return try std.mem.concat(allocator, u8, &[4][]u8{ &header_bytes, self.filename, extra_bytes[0..extra_data_bytes_size], extra_length_bytes[0..extra_length_bytes_size] });
    }
};

pub const GenetaredData = struct {
    prefix: []u8,
    suffix: []u8,
    num_zeroes: usize,

    pub fn deinit(self: *GenetaredData, allocator: std.mem.Allocator) void {
        allocator.free(self.prefix);
        allocator.free(self.suffix);
    }
};

pub const RecordData = union(enum) {
    generated: GenetaredData,
    real_data: []u8,

    pub fn size(self: RecordData) usize {
        return switch (self) {
            RecordData.real_data => |data| data.len,
            RecordData.generated => |data| data.prefix.len + data.suffix.len + data.num_zeroes,
        };
    }

    pub fn deinit(self: *RecordData, allocator: std.mem.Allocator) void {
        switch (self.*) {
            RecordData.real_data => |data| {
                allocator.free(data);
            },
            // RecordData.generated => |*data| data.deinit(allocator),
            RecordData.generated => {},
        }
    }
};

pub const FileRecord = struct {
    header: LocalFileHeader,
    data: RecordData,
    crc_matrix: crc_u.crc_data_type,

    pub fn deinit(self: *FileRecord, allocator: std.mem.Allocator) void {
        self.header.deinit(allocator);
        self.data.deinit(allocator);
    }
};

pub const CompressData = struct {
    data: GenetaredData,
    size: usize,
    crc_matrix: crc_u.crc_data_type,
};

pub const CentralDirectoryHeader = struct {
    offset: usize,
    compressed_size: u64,
    uncompressed_size: u64,
    crc: u32,
    filename: []u8,
    compression_method: u16,

    pub fn new(offset: usize, template: LocalFileHeader) @This() {
        return .{
            .offset = offset,
            .compressed_size = template.compressed_size,
            .uncompressed_size = template.uncompressed_size,
            .crc = template.crc,
            .filename = template.filename,
            .compression_method = template.compression_method,
        };
    }

    pub fn serialize(self: *@This(), zip64: bool, allocator: std.mem.Allocator) ![]u8 {
        const extra_data_size = @bitSizeOf(ExtraDataMin) / 8;

        var extra: ?ExtraDataMin = null;
        if (zip64) {
            extra = ExtraDataMin{
                .tag = 0x0001,
                .size = 8,
                .uncompressed_size = self.uncompressed_size,
            };
        }
        const header = CentralDirectoryData{
            .signature = 0x02014b50,
            .version_by = (0 << 8) | zip_version_num(self.compression_method, zip64),
            .version = zip_version_num(self.compression_method, zip64),
            .flags = 0,
            .compression_method = self.compression_method,
            .mod_time = mod_time,
            .mod_date = mod_date,
            .crc = self.crc,
            .compressed_size = @intCast(self.compressed_size),
            .uncompressed_size = if (zip64) @as(u32, 0xffffffff) else @as(u32, @intCast(self.uncompressed_size)),
            .filename_length = @intCast(self.filename.len),
            .extra_length = if (extra != null) @as(u16, extra_data_size) else 0,
            .comment_size = 0,
            .start_disk = 0,
            .internal_attributes = 0,
            .external_attributes = 0,
            .offset = @intCast(self.offset),
        };
        var header_bytes: [@bitSizeOf(CentralDirectoryData) / 8]u8 = @bitCast(header);
        var extra_bytes: [extra_data_size]u8 = if (extra != null) @bitCast(extra.?) else [_]u8{0} ** (extra_data_size);
        const extra_data_bytes_size: usize = if (extra != null) extra_data_size else 0;
        if (comptime native_endian != .little) {
            std.mem.reverse(u8, &header_bytes);

            if (zip64) {
                std.mem.reverse(u8, &extra_bytes);
            }
        }
        return try std.mem.concat(allocator, u8, &[3][]u8{ &header_bytes, self.filename, extra_bytes[0..extra_data_bytes_size] });
    }
};

const Zip64EOCDRecord = packed struct {
    signature: u32,
    remainder: u64,
    version_by: u16,
    version: u16,
    num_disk: u32,
    cd_disk: u32,
    num_entries_disk: u64,
    num_entries: u64,
    cd_size: u64,
    cd_offset: u64,

    pub fn new(signature: u32, remainder: u64, version_by: u16, version: u16, num_disk: u32, cd_disk: u32, num_entries_disk: u64, num_entries: u64, cd_size: u64, cd_offset: u64) @This() {
        return .{
            .cd_offset = cd_offset,
            .cd_size = cd_size,
            .num_entries = num_entries,
            .num_entries_disk = num_entries_disk,
            .cd_disk = cd_disk,
            .num_disk = num_disk,
            .version = version,
            .version_by = version_by,
            .remainder = remainder,
            .signature = signature,
        };
    }
};

const Zip64EOCDLocator = packed struct {
    signature: u32,
    num_disk_zip64: u32,
    offset: u64,
    total_disk: u32,

    pub fn new(signature: u32, num_disk_zip64: u32, offset: u64, total_disk: u32) @This() {
        return .{
            .total_disk = total_disk,
            .offset = offset,
            .num_disk_zip64 = num_disk_zip64,
            .signature = signature,
        };
    }
};

const EOCDRecord = packed struct {
    signature: u32,
    num_disk: u16,
    cd_disk: u16,
    cd_num_entries: u16,
    cd_num_entries_total: u16,
    cd_size: u32,
    cd_offset: u32,
    comment_length: u16,

    pub fn new(signature: u32, num_disk: u16, cd_disk: u16, cd_num_entries: u16, cd_num_entries_total: u16, cd_size: u32, cd_offset: u32, comment_length: u16) @This() {
        return .{
            .comment_length = comment_length,
            .cd_offset = cd_offset,
            .cd_size = cd_size,
            .cd_num_entries_total = cd_num_entries_total,
            .cd_num_entries = cd_num_entries,
            .cd_disk = cd_disk,
            .num_disk = num_disk,
            .signature = signature,
        };
    }
};

pub const EndOfCentralDirectory = struct {
    num_entries: usize,
    cd_size: usize,
    cd_offset: usize,

    pub fn new(num_entries: usize, cd_size: usize, cd_offset: usize) @This() {
        return .{
            .num_entries = num_entries,
            .cd_size = cd_size,
            .cd_offset = cd_offset,
        };
    }

    pub fn serialize(self: *@This(), zip64: bool, method: u16, allocator: std.mem.Allocator) ![]u8 {
        const zip64_eocd_record_size = @bitSizeOf(Zip64EOCDRecord) / 8;
        const zip64_eocd_locator_size = @bitSizeOf(Zip64EOCDLocator) / 8;
        const eocd_record_size = @bitSizeOf(EOCDRecord) / 8;

        var zip64_eocd_record: ?Zip64EOCDRecord = null;
        var zip64_eocd_locator: ?Zip64EOCDLocator = null;

        if (zip64) {
            zip64_eocd_record = Zip64EOCDRecord.new(
                0x06064b50,
                44,
                zip_version_num(method, zip64),
                zip_version_num(method, zip64),
                0,
                0,
                self.num_entries,
                self.num_entries,
                self.cd_size,
                self.cd_offset,
            );

            zip64_eocd_locator = Zip64EOCDLocator.new(
                0x07064b50,
                0,
                self.cd_offset + self.cd_size,
                1,
            );
        }

        const eocd_record = EOCDRecord.new(
            0x06054b50,
            0,
            0,
            if (zip64) 0xffff else @as(u16, @intCast(self.num_entries)),
            if (zip64) 0xffff else @as(u16, @intCast(self.num_entries)),
            if (zip64) 0xffffffff else @as(u32, @intCast(self.cd_size)),
            if (zip64) 0xffffffff else @as(u32, @intCast(self.cd_offset)),
            0,
        );
        var eocd_record_bytes: [eocd_record_size]u8 = @bitCast(eocd_record);
        var zip64_eocd_record_bytes: [zip64_eocd_record_size]u8 = if (zip64_eocd_record != null) @bitCast(zip64_eocd_record.?) else [_]u8{0} ** (zip64_eocd_record_size);
        var zip64_eocd_locator_bytes: [zip64_eocd_locator_size]u8 = if (zip64_eocd_locator != null) @bitCast(zip64_eocd_locator.?) else [_]u8{0} ** (zip64_eocd_locator_size);
        const trim_zip64_eocd_record_size: usize = if (zip64_eocd_record != null) zip64_eocd_record_size else 0;
        const trim_zip64_eocd_locator_size: usize = if (zip64_eocd_locator != null) zip64_eocd_locator_size else 0;

        if (comptime native_endian != .little) {
            std.mem.reverse(u8, &eocd_record_bytes);

            if (zip64) {
                std.mem.reverse(u8, &zip64_eocd_record_bytes);
                std.mem.reverse(u8, &zip64_eocd_locator_bytes);
            }
        }
        return try std.mem.concat(allocator, u8, &[3][]u8{ zip64_eocd_record_bytes[0..trim_zip64_eocd_record_size], zip64_eocd_locator_bytes[0..trim_zip64_eocd_locator_size], &eocd_record_bytes });
    }
};
