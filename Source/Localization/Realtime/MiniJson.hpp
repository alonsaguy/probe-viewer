#pragma once

// A minimal JSON reader (objects, arrays, numbers, strings, booleans, null) for the localizer
// bundle's manifest. No dependency beyond the standard library, so the localizer core can be
// built and tested without JUCE.

#include <cctype>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
class Json
{
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string text;
    std::vector<Json> items;
    std::map<std::string, Json> members;

    static Json parse (const std::string& s)
    {
        size_t i = 0;
        Json v = parseValue (s, i);
        skip (s, i);
        if (i != s.size())
            throw std::runtime_error ("JSON: trailing characters");
        return v;
    }

    bool has (const std::string& key) const { return type == Type::Object && members.count (key) > 0; }

    const Json& operator[] (const std::string& key) const
    {
        auto it = members.find (key);
        if (type != Type::Object || it == members.end())
            throw std::runtime_error ("JSON: missing key '" + key + "'");
        return it->second;
    }

    const Json& operator[] (size_t index) const
    {
        if (type != Type::Array || index >= items.size())
            throw std::runtime_error ("JSON: index out of range");
        return items[index];
    }

    size_t size() const { return type == Type::Array ? items.size() : members.size(); }

    double num() const
    {
        if (type == Type::Bool)
            return boolean ? 1.0 : 0.0;
        if (type != Type::Number)
            throw std::runtime_error ("JSON: not a number");
        return number;
    }

    const std::string& str() const
    {
        if (type != Type::String)
            throw std::runtime_error ("JSON: not a string");
        return text;
    }

    std::vector<double> numbers() const
    {
        std::vector<double> out;
        for (const auto& item : items)
            out.push_back (item.num());
        return out;
    }

private:
    static void skip (const std::string& s, size_t& i)
    {
        while (i < s.size() && std::isspace (static_cast<unsigned char> (s[i])))
            ++i;
    }

    static Json parseValue (const std::string& s, size_t& i)
    {
        skip (s, i);
        if (i >= s.size())
            throw std::runtime_error ("JSON: unexpected end");
        const char c = s[i];
        Json v;
        if (c == '{')
        {
            v.type = Type::Object;
            ++i;
            skip (s, i);
            if (s[i] == '}')
            {
                ++i;
                return v;
            }
            while (true)
            {
                skip (s, i);
                std::string key = parseString (s, i);
                skip (s, i);
                if (s[i++] != ':')
                    throw std::runtime_error ("JSON: expected ':'");
                v.members[key] = parseValue (s, i);
                skip (s, i);
                if (s[i] == ',')
                {
                    ++i;
                    continue;
                }
                if (s[i++] != '}')
                    throw std::runtime_error ("JSON: expected '}'");
                return v;
            }
        }
        if (c == '[')
        {
            v.type = Type::Array;
            ++i;
            skip (s, i);
            if (s[i] == ']')
            {
                ++i;
                return v;
            }
            while (true)
            {
                v.items.push_back (parseValue (s, i));
                skip (s, i);
                if (s[i] == ',')
                {
                    ++i;
                    continue;
                }
                if (s[i++] != ']')
                    throw std::runtime_error ("JSON: expected ']'");
                return v;
            }
        }
        if (c == '"')
        {
            v.type = Type::String;
            v.text = parseString (s, i);
            return v;
        }
        if (s.compare (i, 4, "true") == 0)
        {
            v.type = Type::Bool;
            v.boolean = true;
            i += 4;
            return v;
        }
        if (s.compare (i, 5, "false") == 0)
        {
            v.type = Type::Bool;
            i += 5;
            return v;
        }
        if (s.compare (i, 4, "null") == 0)
        {
            i += 4;
            return v;
        }
        if (s.compare (i, 3, "NaN") == 0)
        {
            v.type = Type::Number;
            v.number = std::numeric_limits<double>::quiet_NaN();
            i += 3;
            return v;
        }
        // number
        const char* begin = s.c_str() + i;
        char* end = nullptr;
        v.type = Type::Number;
        v.number = std::strtod (begin, &end);
        if (end == begin)
            throw std::runtime_error ("JSON: bad value");
        i += static_cast<size_t> (end - begin);
        return v;
    }

    static std::string parseString (const std::string& s, size_t& i)
    {
        if (s[i] != '"')
            throw std::runtime_error ("JSON: expected string");
        ++i;
        std::string out;
        while (i < s.size() && s[i] != '"')
        {
            if (s[i] == '\\' && i + 1 < s.size())
            {
                const char e = s[i + 1];
                i += 2;
                switch (e)
                {
                    case 'n': out += '\n'; break;
                    case 't': out += '\t'; break;
                    case 'r': out += '\r'; break;
                    case 'b': out += '\b'; break;
                    case 'f': out += '\f'; break;
                    case 'u':
                    {
                        // keep ASCII; replace other code points by '?'
                        const unsigned code = static_cast<unsigned> (std::strtoul (s.substr (i, 4).c_str(), nullptr, 16));
                        out += code < 128 ? static_cast<char> (code) : '?';
                        i += 4;
                        break;
                    }
                    default: out += e; break;
                }
                continue;
            }
            out += s[i++];
        }
        ++i;  // closing quote
        return out;
    }
};
} // namespace ProbeViewer::Localization::Realtime
