#include "warp_hypad/core/json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace warp_hypad {
namespace {

// --- Writing ---------------------------------------------------------------

// Shortest representation that reparses to the identical double.
std::string number_to_string(double v) {
    if (!std::isfinite(v)) {
        // JSON has no NaN/Infinity. Emitting null keeps the document parseable
        // and makes the loss loud rather than silent.
        return "null";
    }
    char buf[64];
    for (int prec = 15; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, v);
        if (std::strtod(buf, nullptr) == v) return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

void escape_string(const std::string& s, std::ostringstream& os) {
    os << '"';
    for (char c : s) {
        switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\b': os << "\\b";  break;
            case '\f': os << "\\f";  break;
            case '\n': os << "\\n";  break;
            case '\r': os << "\\r";  break;
            case '\t': os << "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  static_cast<unsigned char>(c));
                    os << buf;
                } else {
                    os << c;
                }
        }
    }
    os << '"';
}

void dump_value(const JsonValue& v, int indent, int depth,
                std::ostringstream& os) {
    const bool pretty = indent > 0;
    const std::string pad(pretty ? static_cast<std::size_t>(indent * (depth + 1)) : 0, ' ');
    const std::string pad_close(pretty ? static_cast<std::size_t>(indent * depth) : 0, ' ');
    const char* nl = pretty ? "\n" : "";
    const char* sp = pretty ? " " : "";

    switch (v.type()) {
        case JsonValue::Type::Null:   os << "null"; break;
        case JsonValue::Type::Bool:   os << (v.as_bool() ? "true" : "false"); break;
        case JsonValue::Type::Number: os << number_to_string(v.as_number()); break;
        case JsonValue::Type::String: escape_string(v.as_string(), os); break;
        case JsonValue::Type::Array: {
            if (v.items().empty()) { os << "[]"; break; }
            os << '[' << nl;
            for (std::size_t i = 0; i < v.items().size(); ++i) {
                os << pad;
                dump_value(v.items()[i], indent, depth + 1, os);
                if (i + 1 < v.items().size()) os << ',';
                os << nl;
            }
            os << pad_close << ']';
            break;
        }
        case JsonValue::Type::Object: {
            if (v.members().empty()) { os << "{}"; break; }
            os << '{' << nl;
            std::size_t i = 0;
            for (const auto& kv : v.members()) {
                os << pad;
                escape_string(kv.first, os);
                os << ':' << sp;
                dump_value(kv.second, indent, depth + 1, os);
                if (++i < v.members().size()) os << ',';
                os << nl;
            }
            os << pad_close << '}';
            break;
        }
    }
}

// --- Parsing ---------------------------------------------------------------

class Parser {
public:
    Parser(const std::string& text) : s_(text) {}

    bool parse(JsonValue& out, std::string& error) {
        skip_ws();
        if (!parse_value(out)) { error = error_; return false; }
        skip_ws();
        if (pos_ != s_.size()) {
            return fail("trailing content after top-level value", error);
        }
        return true;
    }

private:
    const std::string& s_;
    std::size_t pos_ = 0;
    std::string error_;

    bool fail(const std::string& msg) {
        if (error_.empty()) {
            std::ostringstream os;
            os << msg << " at byte offset " << pos_;
            error_ = os.str();
        }
        return false;
    }
    bool fail(const std::string& msg, std::string& error) {
        fail(msg);
        error = error_;
        return false;
    }

    bool eof() const { return pos_ >= s_.size(); }
    char peek() const { return s_[pos_]; }

    void skip_ws() {
        while (!eof()) {
            const char c = s_[pos_];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') ++pos_;
            else break;
        }
    }

    bool literal(const char* lit) {
        const std::size_t n = std::char_traits<char>::length(lit);
        if (s_.compare(pos_, n, lit) != 0) return false;
        pos_ += n;
        return true;
    }

    bool parse_value(JsonValue& out) {
        if (eof()) return fail("unexpected end of input");
        switch (peek()) {
            case 'n':
                if (!literal("null")) return fail("invalid literal");
                out = JsonValue();
                return true;
            case 't':
                if (!literal("true")) return fail("invalid literal");
                out = JsonValue(true);
                return true;
            case 'f':
                if (!literal("false")) return fail("invalid literal");
                out = JsonValue(false);
                return true;
            case '"': {
                std::string str;
                if (!parse_string(str)) return false;
                out = JsonValue(std::move(str));
                return true;
            }
            case '[': return parse_array(out);
            case '{': return parse_object(out);
            default:  return parse_number(out);
        }
    }

    bool parse_string(std::string& out) {
        if (eof() || peek() != '"') return fail("expected '\"'");
        ++pos_;
        out.clear();
        while (true) {
            if (eof()) return fail("unterminated string");
            const char c = s_[pos_++];
            if (c == '"') return true;
            if (c != '\\') { out.push_back(c); continue; }
            if (eof()) return fail("unterminated escape");
            const char e = s_[pos_++];
            switch (e) {
                case '"':  out.push_back('"');  break;
                case '\\': out.push_back('\\'); break;
                case '/':  out.push_back('/');  break;
                case 'b':  out.push_back('\b'); break;
                case 'f':  out.push_back('\f'); break;
                case 'n':  out.push_back('\n'); break;
                case 'r':  out.push_back('\r'); break;
                case 't':  out.push_back('\t'); break;
                case 'u': {
                    if (pos_ + 4 > s_.size()) return fail("truncated \\u escape");
                    unsigned code = 0;
                    for (int i = 0; i < 4; ++i) {
                        const char h = s_[pos_ + static_cast<std::size_t>(i)];
                        unsigned d;
                        if (h >= '0' && h <= '9') d = static_cast<unsigned>(h - '0');
                        else if (h >= 'a' && h <= 'f') d = static_cast<unsigned>(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F') d = static_cast<unsigned>(h - 'A' + 10);
                        else return fail("bad hex digit in \\u escape");
                        code = code * 16 + d;
                    }
                    pos_ += 4;
                    // UTF-8 encode. Surrogate pairs are not reassembled; the
                    // project writes ASCII keys and values, and mangling a
                    // surrogate silently would be worse than not handling it.
                    if (code < 0x80) {
                        out.push_back(static_cast<char>(code));
                    } else if (code < 0x800) {
                        out.push_back(static_cast<char>(0xC0 | (code >> 6)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    } else {
                        out.push_back(static_cast<char>(0xE0 | (code >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    break;
                }
                default: return fail("unknown escape character");
            }
        }
    }

    bool parse_number(JsonValue& out) {
        const std::size_t start = pos_;
        if (!eof() && (peek() == '-' || peek() == '+')) ++pos_;
        bool any_digit = false;
        while (!eof() && peek() >= '0' && peek() <= '9') { ++pos_; any_digit = true; }
        if (!eof() && peek() == '.') {
            ++pos_;
            while (!eof() && peek() >= '0' && peek() <= '9') { ++pos_; any_digit = true; }
        }
        if (!any_digit) return fail("expected a value");
        if (!eof() && (peek() == 'e' || peek() == 'E')) {
            ++pos_;
            if (!eof() && (peek() == '-' || peek() == '+')) ++pos_;
            bool exp_digit = false;
            while (!eof() && peek() >= '0' && peek() <= '9') { ++pos_; exp_digit = true; }
            if (!exp_digit) return fail("malformed exponent");
        }
        const std::string tok = s_.substr(start, pos_ - start);
        out = JsonValue(std::strtod(tok.c_str(), nullptr));
        return true;
    }

    bool parse_array(JsonValue& out) {
        ++pos_;  // '['
        out = JsonValue::array();
        skip_ws();
        if (!eof() && peek() == ']') { ++pos_; return true; }
        while (true) {
            skip_ws();
            JsonValue item;
            if (!parse_value(item)) return false;
            out.push_back(std::move(item));
            skip_ws();
            if (eof()) return fail("unterminated array");
            if (peek() == ',') { ++pos_; continue; }
            if (peek() == ']') { ++pos_; return true; }
            return fail("expected ',' or ']'");
        }
    }

    bool parse_object(JsonValue& out) {
        ++pos_;  // '{'
        out = JsonValue::object();
        skip_ws();
        if (!eof() && peek() == '}') { ++pos_; return true; }
        while (true) {
            skip_ws();
            std::string key;
            if (!parse_string(key)) return false;
            skip_ws();
            if (eof() || peek() != ':') return fail("expected ':'");
            ++pos_;
            skip_ws();
            JsonValue val;
            if (!parse_value(val)) return false;
            out.set(key, std::move(val));
            skip_ws();
            if (eof()) return fail("unterminated object");
            if (peek() == ',') { ++pos_; continue; }
            if (peek() == '}') { ++pos_; return true; }
            return fail("expected ',' or '}'");
        }
    }
};

}  // namespace

bool json_parse(const std::string& text, JsonValue& out, std::string& error) {
    Parser p(text);
    return p.parse(out, error);
}

std::string json_dump(const JsonValue& v, int indent) {
    std::ostringstream os;
    dump_value(v, indent, 0, os);
    if (indent > 0) os << '\n';
    return os.str();
}

bool json_read_file(const std::string& path, JsonValue& out,
                    std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        error = "cannot open file for reading: " + path;
        return false;
    }
    std::ostringstream buf;
    buf << in.rdbuf();
    return json_parse(buf.str(), out, error);
}

bool json_write_file(const std::string& path, const JsonValue& v,
                     std::string& error, int indent) {
    std::ofstream os(path, std::ios::binary);
    if (!os) {
        error = "cannot open file for writing: " + path;
        return false;
    }
    os << json_dump(v, indent);
    if (!os) {
        error = "write failed: " + path;
        return false;
    }
    return true;
}

}  // namespace warp_hypad
