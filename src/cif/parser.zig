//! CIF parser: reads a CIF source string into a Document.

const std = @import("std");
const Allocator = std.mem.Allocator;

const tokenizer_mod = @import("tokenizer.zig");
const Tokenizer = tokenizer_mod.Tokenizer;
const Token = tokenizer_mod.Token;
const TokenType = tokenizer_mod.TokenType;

const types = @import("types.zig");
pub const Document = types.Document;
pub const Block = types.Block;
pub const Loop = types.Loop;
pub const Item = types.Item;
pub const Pair = types.Pair;

/// Parse a CIF source string and return a Document.
/// The caller owns the returned Document and must call deinit().
/// All string slices in the Document point into `source`.
pub fn readString(allocator: Allocator, source: []const u8) !Document {
    var doc = Document.init(allocator);
    errdefer doc.deinit();

    var tok = try Tokenizer.init(source);

    // Pending token: a token consumed during loop value reading that was not
    // a value and needs to be re-processed in the main loop.
    var pending: ?Token = null;

    while (true) {
        const token = if (pending) |p| blk: {
            pending = null;
            break :blk p;
        } else tok.next();

        switch (token.type) {
            .eof => break,

            .data => {
                // Block name is everything after "data_"
                const block_name = token.text(source)[5..];
                try doc.blocks.append(allocator, Block{ .name = block_name });
            },

            .tag => {
                if (doc.blocks.items.len == 0) return error.TagOutsideBlock;
                const block = &doc.blocks.items[doc.blocks.items.len - 1];
                const tag_text = token.text(source);

                // The next token must be a value
                const val_tok = tok.next();
                if (val_tok.type != .value) return error.ExpectedValue;
                const val_text = val_tok.text(source);

                try block.items.append(allocator, .{ .pair = Pair{
                    .tag = tag_text,
                    .value = val_text,
                } });
            },

            .loop => {
                if (doc.blocks.items.len == 0) return error.LoopOutsideBlock;
                const block = &doc.blocks.items[doc.blocks.items.len - 1];

                var loop = Loop{};
                errdefer loop.deinit(allocator);

                // Read tags
                var first_non_tag: Token = undefined;
                while (true) {
                    const t = tok.next();
                    if (t.type == .tag) {
                        try loop.tags.append(allocator, t.text(source));
                    } else {
                        first_non_tag = t;
                        break;
                    }
                }

                // Read values — the first non-tag token may already be a value
                var next_tok = first_non_tag;
                while (next_tok.type == .value) {
                    try loop.values.append(allocator, next_tok.text(source));
                    next_tok = tok.next();
                }
                // next_tok is the first non-value token after the loop body
                // Store it as pending so the main loop processes it
                if (next_tok.type != .eof) {
                    pending = next_tok;
                }

                // Validate that value count is divisible by tag count.
                if (loop.tags.items.len > 0 and
                    loop.values.items.len % loop.tags.items.len != 0)
                {
                    return error.LoopValueCountMismatch;
                }

                try block.items.append(allocator, .{ .loop = loop });
            },

            // save frames: skip for now (CIF dictionaries use these)
            .save_begin, .save_end => {},

            // Unterminated quoted string or text field
            .invalid => return error.UnterminatedValue,

            .value => {
                // Bare value outside a loop/pair context indicates malformed CIF
                return error.UnexpectedValue;
            },
        }
    }

    return doc;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "parse simple document" {
    const source =
        \\data_TEST
        \\_entry.id TEST
        \\loop_
        \\_atom.x
        \\_atom.y
        \\1.0 2.0
        \\3.0 4.0
    ;

    var doc = try readString(testing.allocator, source);
    defer doc.deinit();

    // One block named "TEST"
    try testing.expectEqual(@as(usize, 1), doc.blocks.items.len);

    const block = doc.findBlock("TEST") orelse return error.BlockNotFound;

    // findValue
    const entry_id = block.findValue("_entry.id") orelse return error.ValueNotFound;
    try testing.expectEqualStrings("TEST", entry_id);

    // findLoop
    const loop = block.findLoop("_atom.x") orelse return error.LoopNotFound;
    try testing.expectEqual(@as(usize, 2), loop.width());
    try testing.expectEqual(@as(usize, 2), loop.length());

    const v00 = loop.val(0, 0) orelse return error.ValNotFound;
    try testing.expectEqualStrings("1.0", v00);

    const v11 = loop.val(1, 1) orelse return error.ValNotFound;
    try testing.expectEqualStrings("4.0", v11);
}

test "parse multiple blocks" {
    const source =
        \\data_FIRST
        \\_key value1
        \\data_SECOND
        \\_key value2
    ;

    var doc = try readString(testing.allocator, source);
    defer doc.deinit();

    try testing.expectEqual(@as(usize, 2), doc.blocks.items.len);

    const b1 = doc.findBlock("FIRST") orelse return error.BlockNotFound;
    try testing.expectEqualStrings("value1", b1.findValue("_key").?);

    const b2 = doc.findBlock("SECOND") orelse return error.BlockNotFound;
    try testing.expectEqualStrings("value2", b2.findValue("_key").?);
}

test "parse loop followed by pair" {
    const source =
        \\data_X
        \\loop_
        \\_a
        \\1 2 3
        \\_after loop_end
    ;

    var doc = try readString(testing.allocator, source);
    defer doc.deinit();

    const block = doc.findBlock("X") orelse return error.BlockNotFound;
    const loop = block.findLoop("_a") orelse return error.LoopNotFound;
    try testing.expectEqual(@as(usize, 3), loop.length());

    const after = block.findValue("_after") orelse return error.ValueNotFound;
    try testing.expectEqualStrings("loop_end", after);
}
