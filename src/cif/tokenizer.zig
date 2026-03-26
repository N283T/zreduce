//! CIF tokenizer: zero-copy token stream (tokens are slices into source).

const char_table = @import("char_table.zig");

pub const TokenType = enum {
    data, // data_BLOCKNAME
    loop, // loop_
    tag, // _tag_name
    value, // unquoted, single-quoted, double-quoted, or semicolon-delimited
    save_begin, // save_FRAMENAME
    save_end, // save_
    invalid, // unterminated quoted string or text field
    eof,
};

pub const Token = struct {
    type: TokenType,
    start: u32,
    end: u32,

    pub fn text(self: Token, source: []const u8) []const u8 {
        return source[self.start..self.end];
    }
};

pub const Tokenizer = struct {
    source: []const u8,
    pos: u32,

    pub fn init(source: []const u8) Tokenizer {
        return .{ .source = source, .pos = 0 };
    }

    pub fn next(self: *Tokenizer) Token {
        self.skipWhitespaceAndComments();

        if (self.pos >= self.source.len) {
            return .{ .type = .eof, .start = self.pos, .end = self.pos };
        }

        const c = self.source[self.pos];

        // Semicolon text field: only when at start of line
        if (c == ';' and self.isStartOfLine()) {
            return self.readTextField();
        }

        // Quoted strings
        if (c == '\'' or c == '"') {
            return self.readQuotedString(c);
        }

        // Tags
        if (c == '_') {
            return self.readTag();
        }

        // Bare word: loop_, data_*, save_*, or plain value
        const start = self.pos;
        self.advanceToWhitespace();
        const end = self.pos;
        const word = self.source[start..end];

        if (std.ascii.eqlIgnoreCase(word, "loop_")) {
            return .{ .type = .loop, .start = start, .end = end };
        }

        if (word.len >= 5 and std.ascii.startsWithIgnoreCase(word, "data_")) {
            return .{ .type = .data, .start = start, .end = end };
        }

        if (word.len >= 5 and std.ascii.startsWithIgnoreCase(word, "save_")) {
            // Bare "save_" (length 5) marks end of save frame; longer is save_begin
            if (word.len == 5) {
                return .{ .type = .save_end, .start = start, .end = end };
            }
            return .{ .type = .save_begin, .start = start, .end = end };
        }

        return .{ .type = .value, .start = start, .end = end };
    }

    fn skipWhitespaceAndComments(self: *Tokenizer) void {
        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (char_table.isWhitespace(c)) {
                self.pos += 1;
            } else if (c == '#') {
                // Skip to end of line
                while (self.pos < self.source.len and self.source[self.pos] != '\n') {
                    self.pos += 1;
                }
            } else {
                break;
            }
        }
    }

    fn isStartOfLine(self: *const Tokenizer) bool {
        if (self.pos == 0) return true;
        return self.source[self.pos - 1] == '\n';
    }

    /// Read a semicolon-delimited text field.
    /// The opening ';' must be at the start of a line.
    /// Content runs until a bare ';' at the start of a new line.
    /// The returned token covers just the inner content (after opening ';' newline,
    /// before closing ';'), consistent with how CIF text fields are used.
    fn readTextField(self: *Tokenizer) Token {
        // Skip the opening ';'
        self.pos += 1;
        // Skip the newline immediately after the opening semicolon (if present)
        if (self.pos < self.source.len and self.source[self.pos] == '\n') {
            self.pos += 1;
        } else if (self.pos < self.source.len and self.source[self.pos] == '\r') {
            self.pos += 1;
            if (self.pos < self.source.len and self.source[self.pos] == '\n') {
                self.pos += 1;
            }
        }

        const content_start = self.pos;

        // Scan for closing ';' at start of line
        while (self.pos < self.source.len) {
            if (self.source[self.pos] == ';' and self.isStartOfLine()) {
                const content_end = self.pos;
                // Skip the closing ';'
                self.pos += 1;
                return .{ .type = .value, .start = content_start, .end = content_end };
            }
            self.pos += 1;
        }

        // Unterminated text field
        return .{ .type = .invalid, .start = content_start, .end = self.pos };
    }

    /// Read a single- or double-quoted string.
    /// The token covers only the inner content (excluding quotes).
    fn readQuotedString(self: *Tokenizer, quote: u8) Token {
        // Skip opening quote
        self.pos += 1;
        const content_start = self.pos;

        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (c == quote) {
                // Quote must be followed by whitespace or EOF to close the string
                const next_pos = self.pos + 1;
                const at_end = next_pos >= self.source.len;
                const next_is_ws = !at_end and char_table.isWhitespace(self.source[next_pos]);
                if (at_end or next_is_ws) {
                    const content_end = self.pos;
                    self.pos += 1; // Skip closing quote
                    return .{ .type = .value, .start = content_start, .end = content_end };
                }
            }
            self.pos += 1;
        }

        // Unterminated quoted string
        return .{ .type = .invalid, .start = content_start, .end = self.pos };
    }

    fn readTag(self: *Tokenizer) Token {
        const start = self.pos;
        self.advanceToWhitespace();
        return .{ .type = .tag, .start = start, .end = self.pos };
    }

    fn advanceToWhitespace(self: *Tokenizer) void {
        while (self.pos < self.source.len and !char_table.isWhitespace(self.source[self.pos])) {
            self.pos += 1;
        }
    }
};

// --- Tests ---

const std = @import("std");

test "tokenize data block" {
    const source = "data_TEST\n_tag value\n";
    var tok = Tokenizer.init(source);
    const first = tok.next();
    try std.testing.expectEqual(TokenType.data, first.type);
    try std.testing.expectEqualStrings("data_TEST", first.text(source));
}

test "tokenize loop" {
    const source = "loop_\n_col1\n_col2\nA B\n";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.loop, t0.type);

    const t1 = tok.next();
    try std.testing.expectEqual(TokenType.tag, t1.type);

    const t2 = tok.next();
    try std.testing.expectEqual(TokenType.tag, t2.type);

    const t3 = tok.next();
    try std.testing.expectEqual(TokenType.value, t3.type);

    const t4 = tok.next();
    try std.testing.expectEqual(TokenType.value, t4.type);

    const t5 = tok.next();
    try std.testing.expectEqual(TokenType.eof, t5.type);
}

test "tokenize quoted strings" {
    const source = "'hello world' \"another string\"";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.value, t0.type);
    try std.testing.expectEqualStrings("hello world", t0.text(source));

    const t1 = tok.next();
    try std.testing.expectEqual(TokenType.value, t1.type);
    try std.testing.expectEqualStrings("another string", t1.text(source));

    const t2 = tok.next();
    try std.testing.expectEqual(TokenType.eof, t2.type);
}

test "tokenize semicolon text field" {
    const source = ";\nmulti\nline\nvalue\n;\n";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.value, t0.type);
    try std.testing.expectEqualStrings("multi\nline\nvalue\n", t0.text(source));

    const t1 = tok.next();
    try std.testing.expectEqual(TokenType.eof, t1.type);
}

test "skip comments" {
    const source = "# comment\ndata_X\n# another\n_tag val";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.data, t0.type);
    try std.testing.expectEqualStrings("data_X", t0.text(source));

    const t1 = tok.next();
    try std.testing.expectEqual(TokenType.tag, t1.type);
    try std.testing.expectEqualStrings("_tag", t1.text(source));

    const t2 = tok.next();
    try std.testing.expectEqual(TokenType.value, t2.type);
    try std.testing.expectEqualStrings("val", t2.text(source));

    const t3 = tok.next();
    try std.testing.expectEqual(TokenType.eof, t3.type);
}

test "tokenize save frame" {
    const source = "save_myframe\n_tag val\nsave_\n";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.save_begin, t0.type);
    try std.testing.expectEqualStrings("save_myframe", t0.text(source));

    _ = tok.next(); // _tag
    _ = tok.next(); // val

    const t3 = tok.next();
    try std.testing.expectEqual(TokenType.save_end, t3.type);
    try std.testing.expectEqualStrings("save_", t3.text(source));
}

test "unterminated quoted string returns invalid" {
    const source = "'hello world";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.invalid, t0.type);
}

test "unterminated text field returns invalid" {
    const source = ";\nsome text without closing\n";
    var tok = Tokenizer.init(source);

    const t0 = tok.next();
    try std.testing.expectEqual(TokenType.invalid, t0.type);
}
