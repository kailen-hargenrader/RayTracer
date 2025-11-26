#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <variant>
#include <optional>

namespace rt {

class Json {
public:
    struct Value;
    using Array = std::vector<Value>;
    using Object = std::unordered_map<std::string, Value>;
    using Variant = std::variant<std::nullptr_t, bool, double, std::string, Array, Object>;

    struct Value {
        Variant v;
        Value() : v(nullptr) {}
        Value(std::nullptr_t) : v(nullptr) {}
        Value(bool b) : v(b) {}
        Value(double d) : v(d) {}
        Value(const std::string& s) : v(s) {}
        Value(std::string&& s) : v(std::move(s)) {}
        Value(const Array& a) : v(a) {}
        Value(Array&& a) : v(std::move(a)) {}
        Value(const Object& o) : v(o) {}
        Value(Object&& o) : v(std::move(o)) {}

        bool isNull() const { return std::holds_alternative<std::nullptr_t>(v); }
        bool isBool() const { return std::holds_alternative<bool>(v); }
        bool isNumber() const { return std::holds_alternative<double>(v); }
        bool isString() const { return std::holds_alternative<std::string>(v); }
        bool isArray() const { return std::holds_alternative<Array>(v); }
        bool isObject() const { return std::holds_alternative<Object>(v); }

        const Array& asArray() const { return std::get<Array>(v); }
        const Object& asObject() const { return std::get<Object>(v); }
        const std::string& asString() const { return std::get<std::string>(v); }
        double asNumber() const { return std::get<double>(v); }
        bool asBool() const { return std::get<bool>(v); }

        // Convenience getters with defaults
        double numberOr(double def) const { return isNumber() ? asNumber() : def; }
        std::string stringOr(const std::string& def) const { return isString() ? asString() : def; }
        bool boolOr(bool def) const { return isBool() ? asBool() : def; }

        // Object helpers
        const Value* get(const std::string& key) const {
            if (!isObject()) return nullptr;
            const auto& o = asObject();
            auto it = o.find(key);
            return it == o.end() ? nullptr : &it->second;
        }
    };

    static bool parse(const std::string& s, Value& outValue, std::string* errorMessage = nullptr);
    static bool parseFile(const std::string& path, Value& outValue, std::string* errorMessage = nullptr);
};

} // namespace rt


