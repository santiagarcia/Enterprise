#pragma once

#include <map>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// A small, dependency-free JSON value type with a parser and a writer.
//
// Written rather than pulled in because the project needs exactly three things
// from JSON -- read a config, write run metadata, write a verification summary
// -- and a vendored single-header library would be more code than this, would
// need pinning and licence tracking, and would still not give the one property
// that matters here: numbers that round-trip bit-exactly.
//
// Numbers are emitted at the shortest precision that reproduces the exact
// double on re-parse (try 15, then 16, then 17 significant digits). That makes
// the Milestone 9 "write it, read it back, compare exactly" test meaningful
// rather than tolerance-limited, and keeps configs readable: 0.2 stays "0.2"
// instead of becoming "0.20000000000000001".
//
// Object keys are stored in a std::map, so serialisation order is sorted and
// therefore reproducible across runs and platforms. Reproducible bytes matter
// because metadata.json is hashed and compared between runs.
//
// Not supported, deliberately: comments, trailing commas, NaN/Infinity
// literals. A configuration that needs any of those is a configuration that
// will not round-trip.
// ---------------------------------------------------------------------------

namespace warp_hypad {

class JsonValue {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    JsonValue() : type_(Type::Null) {}
    JsonValue(bool b) : type_(Type::Bool), bool_(b) {}
    JsonValue(double n) : type_(Type::Number), num_(n) {}
    JsonValue(int n) : type_(Type::Number), num_(static_cast<double>(n)) {}
    JsonValue(const char* s) : type_(Type::String), str_(s) {}
    JsonValue(std::string s) : type_(Type::String), str_(std::move(s)) {}

    static JsonValue array() {
        JsonValue v;
        v.type_ = Type::Array;
        return v;
    }
    static JsonValue object() {
        JsonValue v;
        v.type_ = Type::Object;
        return v;
    }

    Type type() const { return type_; }
    bool is_null() const { return type_ == Type::Null; }
    bool is_bool() const { return type_ == Type::Bool; }
    bool is_number() const { return type_ == Type::Number; }
    bool is_string() const { return type_ == Type::String; }
    bool is_array() const { return type_ == Type::Array; }
    bool is_object() const { return type_ == Type::Object; }

    bool as_bool() const { return bool_; }
    double as_number() const { return num_; }
    const std::string& as_string() const { return str_; }

    const std::vector<JsonValue>& items() const { return arr_; }
    std::vector<JsonValue>& items() { return arr_; }

    const std::map<std::string, JsonValue>& members() const { return obj_; }
    std::map<std::string, JsonValue>& members() { return obj_; }

    void push_back(JsonValue v) {
        type_ = Type::Array;
        arr_.push_back(std::move(v));
    }

    // Inserts or overwrites. Turns the value into an object if it was not one.
    void set(const std::string& key, JsonValue v) {
        type_ = Type::Object;
        obj_[key] = std::move(v);
    }

    bool has(const std::string& key) const {
        return type_ == Type::Object && obj_.find(key) != obj_.end();
    }

    // Returns nullptr when absent or when this is not an object.
    const JsonValue* find(const std::string& key) const {
        if (type_ != Type::Object) return nullptr;
        auto it = obj_.find(key);
        return it == obj_.end() ? nullptr : &it->second;
    }

private:
    Type type_;
    bool bool_ = false;
    double num_ = 0.0;
    std::string str_;
    std::vector<JsonValue> arr_;
    std::map<std::string, JsonValue> obj_;
};

// Returns false and fills `error` with a message including the byte offset.
bool json_parse(const std::string& text, JsonValue& out, std::string& error);

// Pretty-printed with `indent` spaces per level; indent = 0 emits compact form.
std::string json_dump(const JsonValue& v, int indent = 2);

// Convenience: read/write a whole file. Both return false on I/O failure.
bool json_read_file(const std::string& path, JsonValue& out, std::string& error);
bool json_write_file(const std::string& path, const JsonValue& v,
                     std::string& error, int indent = 2);

}  // namespace warp_hypad
