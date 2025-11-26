#include "raytracer/utils/Json.h"
#include <fstream>
#include <sstream>
#include <cctype>

namespace rt {

namespace {

struct Lexer {
    const char* p;
    const char* end;
    explicit Lexer(const std::string& s) : p(s.data()), end(s.data()+s.size()) {}
    void skipWs() {
        while (p < end && std::isspace(static_cast<unsigned char>(*p))) ++p;
    }
    bool match(char c) {
        skipWs();
        if (p < end && *p == c) { ++p; return true; }
        return false;
    }
    bool peek(char c) {
        skipWs();
        return (p < end && *p == c);
    }
    bool eof() {
        skipWs();
        return p >= end;
    }
};

bool parseString(Lexer& lx, std::string& out) {
    if (!lx.match('\"')) return false;
    std::string s;
    while (lx.p < lx.end) {
        char c = *lx.p++;
        if (c == '\"') { out = std::move(s); return true; }
        if (c == '\\') {
            if (lx.p >= lx.end) return false;
            char e = *lx.p++;
            switch (e) {
                case '\"': s.push_back('\"'); break;
                case '\\': s.push_back('\\'); break;
                case '/': s.push_back('/'); break;
                case 'b': s.push_back('\b'); break;
                case 'f': s.push_back('\f'); break;
                case 'n': s.push_back('\n'); break;
                case 'r': s.push_back('\r'); break;
                case 't': s.push_back('\t'); break;
                default: return false; // no unicode handling
            }
        } else {
            s.push_back(c);
        }
    }
    return false;
}

bool parseNumber(Lexer& lx, double& out) {
    lx.skipWs();
    const char* start = lx.p;
    if (lx.p < lx.end && (*lx.p == '-' || *lx.p == '+')) ++lx.p;
    bool hasDigits = false;
    while (lx.p < lx.end && std::isdigit(static_cast<unsigned char>(*lx.p))) { ++lx.p; hasDigits = true; }
    if (lx.p < lx.end && *lx.p == '.') {
        ++lx.p;
        while (lx.p < lx.end && std::isdigit(static_cast<unsigned char>(*lx.p))) { ++lx.p; hasDigits = true; }
    }
    if (lx.p < lx.end && (*lx.p == 'e' || *lx.p == 'E')) {
        ++lx.p;
        if (lx.p < lx.end && (*lx.p == '+' || *lx.p == '-')) ++lx.p;
        while (lx.p < lx.end && std::isdigit(static_cast<unsigned char>(*lx.p))) ++lx.p;
    }
    if (!hasDigits) return false;
    try {
        out = std::stod(std::string(start, lx.p));
        return true;
    } catch (...) {
        return false;
    }
}

bool parseValue(Lexer& lx, Json::Value& out);

bool parseArray(Lexer& lx, Json::Array& out) {
    if (!lx.match('[')) return false;
    lx.skipWs();
    if (lx.peek(']')) { lx.match(']'); return true; }
    while (true) {
        Json::Value v;
        if (!parseValue(lx, v)) return false;
        out.emplace_back(std::move(v));
        lx.skipWs();
        if (lx.match(']')) return true;
        if (!lx.match(',')) return false;
    }
}

bool parseObject(Lexer& lx, Json::Object& out) {
    if (!lx.match('{')) return false;
    lx.skipWs();
    if (lx.peek('}')) { lx.match('}'); return true; }
    while (true) {
        std::string key;
        if (!parseString(lx, key)) return false;
        if (!lx.match(':')) return false;
        Json::Value v;
        if (!parseValue(lx, v)) return false;
        out.emplace(std::move(key), std::move(v));
        lx.skipWs();
        if (lx.match('}')) return true;
        if (!lx.match(',')) return false;
    }
}

bool parseLiteral(Lexer& lx, const char* lit) {
    while (*lit) {
        if (lx.p >= lx.end || *lx.p != *lit) return false;
        ++lx.p; ++lit;
    }
    return true;
}

bool parseValue(Lexer& lx, Json::Value& out) {
    lx.skipWs();
    if (lx.p >= lx.end) return false;
    char c = *lx.p;
    if (c == '\"') {
        std::string s;
        if (!parseString(lx, s)) return false;
        out = Json::Value(std::move(s));
        return true;
    }
    if (c == '{') {
        Json::Object o;
        if (!parseObject(lx, o)) return false;
        out = Json::Value(std::move(o));
        return true;
    }
    if (c == '[') {
        Json::Array a;
        if (!parseArray(lx, a)) return false;
        out = Json::Value(std::move(a));
        return true;
    }
    if (c == 't') {
        if (!parseLiteral(lx, "true")) return false;
        out = Json::Value(true);
        return true;
    }
    if (c == 'f') {
        if (!parseLiteral(lx, "false")) return false;
        out = Json::Value(false);
        return true;
    }
    if (c == 'n') {
        if (!parseLiteral(lx, "null")) return false;
        out = Json::Value(nullptr);
        return true;
    }
    double num;
    if (parseNumber(lx, num)) { out = Json::Value(num); return true; }
    return false;
}

} // namespace

bool Json::parse(const std::string& s, Value& outValue, std::string* errorMessage) {
    Lexer lx(s);
    if (!parseValue(lx, outValue)) {
        if (errorMessage) *errorMessage = "Failed to parse JSON value.";
        return false;
    }
    if (!lx.eof()) {
        if (errorMessage) *errorMessage = "Extra data after JSON value.";
        return false;
    }
    return true;
}

bool Json::parseFile(const std::string& path, Value& outValue, std::string* errorMessage) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        if (errorMessage) *errorMessage = "Could not open file: " + path;
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return parse(ss.str(), outValue, errorMessage);
}

} // namespace rt


