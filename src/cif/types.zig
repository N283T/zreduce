//! CIF document types: Document, Block, Item, Loop, Pair.

const std = @import("std");
const Allocator = std.mem.Allocator;

// In Zig 0.15, std.ArrayList is the unmanaged variant (no stored allocator).
// append/deinit all take an explicit Allocator argument.

// ---------------------------------------------------------------------------
// Pair
// ---------------------------------------------------------------------------

pub const Pair = struct {
    tag: []const u8,
    value: []const u8,
};

// ---------------------------------------------------------------------------
// Loop
// ---------------------------------------------------------------------------

pub const Loop = struct {
    tags: std.ArrayList([]const u8) = .empty,
    values: std.ArrayList([]const u8) = .empty,

    pub fn deinit(self: *Loop, allocator: Allocator) void {
        self.tags.deinit(allocator);
        self.values.deinit(allocator);
    }

    /// Number of tag columns.
    pub fn width(self: *const Loop) usize {
        return self.tags.items.len;
    }

    /// Number of data rows.
    pub fn length(self: *const Loop) usize {
        const w = self.width();
        if (w == 0) return 0;
        return @divExact(self.values.items.len, w);
    }

    /// Value at (row, col), bounds-checked.
    pub fn val(self: *const Loop, row: usize, col: usize) ?[]const u8 {
        const w = self.width();
        if (w == 0) return null;
        if (col >= w) return null;
        const idx = row * w + col;
        if (idx >= self.values.items.len) return null;
        return self.values.items[idx];
    }

    /// Find column index of tag (case-insensitive). Returns null if not found.
    pub fn findTag(self: *const Loop, tag: []const u8) ?usize {
        for (self.tags.items, 0..) |t, i| {
            if (std.ascii.eqlIgnoreCase(t, tag)) return i;
        }
        return null;
    }
};

// ---------------------------------------------------------------------------
// Item
// ---------------------------------------------------------------------------

pub const ItemType = enum { pair, loop };

pub const Item = union(ItemType) {
    pair: Pair,
    loop: Loop,

    pub fn deinit(self: *Item, allocator: Allocator) void {
        switch (self.*) {
            .pair => {}, // slices are borrowed from source
            .loop => |*l| l.deinit(allocator),
        }
    }
};

// ---------------------------------------------------------------------------
// Block
// ---------------------------------------------------------------------------

pub const Block = struct {
    name: []const u8,
    items: std.ArrayList(Item) = .empty,

    pub fn deinit(self: *Block, allocator: Allocator) void {
        for (self.items.items) |*item| {
            item.deinit(allocator);
        }
        self.items.deinit(allocator);
    }

    /// Find the first loop that contains the given tag (case-insensitive).
    pub fn findLoop(self: *const Block, tag: []const u8) ?*const Loop {
        for (self.items.items) |*item| {
            switch (item.*) {
                .loop => |*l| {
                    if (l.findTag(tag) != null) return l;
                },
                .pair => {},
            }
        }
        return null;
    }

    /// Find the first loop that contains the given tag (mutable).
    pub fn findLoopMut(self: *Block, tag: []const u8) ?*Loop {
        for (self.items.items) |*item| {
            switch (item.*) {
                .loop => |*l| {
                    if (l.findTag(tag) != null) return l;
                },
                .pair => {},
            }
        }
        return null;
    }

    /// Find value of the first pair with the given tag (case-insensitive).
    pub fn findValue(self: *const Block, tag: []const u8) ?[]const u8 {
        for (self.items.items) |*item| {
            switch (item.*) {
                .pair => |p| {
                    if (std.ascii.eqlIgnoreCase(p.tag, tag)) return p.value;
                },
                .loop => {},
            }
        }
        return null;
    }
};

// ---------------------------------------------------------------------------
// Document
// ---------------------------------------------------------------------------

pub const Document = struct {
    allocator: Allocator,
    blocks: std.ArrayList(Block) = .empty,

    pub fn init(allocator: Allocator) Document {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Document) void {
        for (self.blocks.items) |*block| {
            block.deinit(self.allocator);
        }
        self.blocks.deinit(self.allocator);
    }

    /// Find block by name (case-insensitive).
    pub fn findBlock(self: *const Document, name: []const u8) ?*const Block {
        for (self.blocks.items) |*block| {
            if (std.ascii.eqlIgnoreCase(block.name, name)) return block;
        }
        return null;
    }
};
