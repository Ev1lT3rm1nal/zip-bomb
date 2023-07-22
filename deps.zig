// zig fmt: off
const std = @import("std");
const builtin = @import("builtin");
const string = []const u8;
const ModuleDependency = std.build.ModuleDependency;

pub const cache = ".zigmod\\deps";

pub fn addAllTo(exe: *std.build.LibExeObjStep) void {
    checkMinZig(builtin.zig_version, exe);
    const b = exe.step.owner;
    @setEvalBranchQuota(1_000_000);
    for (packages) |pkg| {
        const moddep = pkg.zp(b);
        exe.addModule(moddep.name, moddep.module);
    }
    var llc = false;
    var vcpkg = false;
    inline for (comptime std.meta.declarations(package_data)) |decl| {
        const pkg = @as(Package, @field(package_data, decl.name));
        for (pkg.system_libs) |item| {
            exe.linkSystemLibrary(item);
            llc = true;
        }
        for (pkg.frameworks) |item| {
            if (!builtin.target.isDarwin()) @panic(b.fmt("a dependency is attempting to link to the framework {s}, which is only possible under Darwin", .{item}));
            exe.linkFramework(item);
            llc = true;
        }
        for (pkg.c_include_dirs) |item| {
            exe.addIncludePath(b.fmt("{s}/{s}", .{ @field(dirs, decl.name), item }));
            llc = true;
        }
        for (pkg.c_source_files) |item| {
            exe.addCSourceFile(b.fmt("{s}/{s}", .{ @field(dirs, decl.name), item }), pkg.c_source_flags);
            llc = true;
        }
        vcpkg = vcpkg or pkg.vcpkg;
    }
    if (llc) exe.linkLibC();
    if (builtin.os.tag == .windows and vcpkg) exe.addVcpkgPaths(.static) catch |err| @panic(@errorName(err));
}

pub const Package = struct {
    directory: string,
    pkg: ?Pkg = null,
    c_include_dirs: []const string = &.{},
    c_source_files: []const string = &.{},
    c_source_flags: []const string = &.{},
    system_libs: []const string = &.{},
    frameworks: []const string = &.{},
    vcpkg: bool = false,
    module: ?ModuleDependency = null,

    pub fn zp(self: *Package, b: *std.build.Builder) ModuleDependency {
        var temp: [100]ModuleDependency = undefined;
        const pkg = self.pkg.?;
        for (pkg.dependencies, 0..) |item, i| {
            temp[i] = item.zp(b);
        }
        if (self.module) |mod| {
            return mod;
        }
        const result = ModuleDependency{
            .name = pkg.name,
            .module = b.createModule(.{
                .source_file = pkg.source,
                .dependencies = b.allocator.dupe(ModuleDependency, temp[0..pkg.dependencies.len]) catch @panic("oom"),
            }),
        };
        self.module = result;
        return result;
    }
};

pub const Pkg = struct {
    name: string,
    source: std.build.FileSource,
    dependencies: []const *Package,
};

fn checkMinZig(current: std.SemanticVersion, exe: *std.build.LibExeObjStep) void {
    const min = std.SemanticVersion.parse("null") catch return;
    if (current.order(min).compare(.lt)) @panic(exe.step.owner.fmt("Your Zig version v{} does not meet the minimum build requirement of v{}", .{current, min}));
}

pub const dirs = struct {
    pub const _root = "";
    pub const _ah8gg7veru0u = cache ++ "/../..";
    pub const _xzxo5rnug8wj = cache ++ "/git/github.com/MasterQ32/zig-args";
};

pub const package_data = struct {
    pub var _ah8gg7veru0u = Package{
        .directory = dirs._ah8gg7veru0u,
    };
    pub var _xzxo5rnug8wj = Package{
        .directory = dirs._xzxo5rnug8wj,
        .pkg = Pkg{ .name = "args", .source = .{ .path = dirs._xzxo5rnug8wj ++ "/args.zig" }, .dependencies = &.{} },
    };
    pub var _root = Package{
        .directory = dirs._root,
    };
};

pub const packages = &[_]*Package{
    &package_data._xzxo5rnug8wj,
};

pub const pkgs = struct {
    pub const args = &package_data._xzxo5rnug8wj;
};

pub const imports = struct {
};
