#pragma once

// Reader for NumPy .npy files (format 1.x / 2.x, little-endian, C order), as written by the
// localizer bundle export (ephysatlas.realtime.bundle).

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ProbeViewer::Localization::Realtime
{
struct NpyArray
{
    std::string dtype; // e.g. "<f8", "<f4", "<i4", "<i8", "<u2", "|u1", "|b1"
    std::vector<size_t> shape;
    std::vector<char> bytes;

    size_t count() const
    {
        size_t n = 1;
        for (auto s : shape)
            n *= s;
        return n;
    }

    // Values converted to T, whatever the stored numeric type.
    template <typename T>
    std::vector<T> as() const
    {
        const size_t n = count();
        std::vector<T> out (n);
        const char kind = dtype.size() > 1 ? dtype[1] : '?';
        const int width = dtype.size() > 2 ? std::stoi (dtype.substr (2)) : 1;
        for (size_t i = 0; i < n; ++i)
        {
            const char* p = bytes.data() + i * static_cast<size_t> (width);
            double v = 0.0;
            if (kind == 'f' && width == 8) { double x; std::memcpy (&x, p, 8); v = x; }
            else if (kind == 'f' && width == 4) { float x; std::memcpy (&x, p, 4); v = x; }
            else if (kind == 'i' && width == 8) { int64_t x; std::memcpy (&x, p, 8); v = static_cast<double> (x); }
            else if (kind == 'i' && width == 4) { int32_t x; std::memcpy (&x, p, 4); v = x; }
            else if (kind == 'i' && width == 2) { int16_t x; std::memcpy (&x, p, 2); v = x; }
            else if (kind == 'u' && width == 2) { uint16_t x; std::memcpy (&x, p, 2); v = x; }
            else if (kind == 'u' && width == 1) { uint8_t x; std::memcpy (&x, p, 1); v = x; }
            else if (kind == 'b' && width == 1) { uint8_t x; std::memcpy (&x, p, 1); v = x ? 1.0 : 0.0; }
            else throw std::runtime_error ("npy: unsupported dtype " + dtype);
            out[i] = static_cast<T> (v);
        }
        return out;
    }
};

inline NpyArray loadNpy (const std::string& path)
{
    std::ifstream f (path, std::ios::binary);
    if (! f)
        throw std::runtime_error ("npy: cannot open " + path);
    char magic[6];
    f.read (magic, 6);
    if (std::memcmp (magic, "\x93NUMPY", 6) != 0)
        throw std::runtime_error ("npy: not a .npy file: " + path);
    unsigned char version[2];
    f.read (reinterpret_cast<char*> (version), 2);
    uint32_t headerLen = 0;
    if (version[0] == 1)
    {
        uint16_t h;
        f.read (reinterpret_cast<char*> (&h), 2);
        headerLen = h;
    }
    else
    {
        f.read (reinterpret_cast<char*> (&headerLen), 4);
    }
    std::string header (headerLen, ' ');
    f.read (header.data(), headerLen);

    NpyArray a;
    auto field = [&header] (const std::string& key) -> std::string {
        const auto k = header.find ("'" + key + "'");
        if (k == std::string::npos)
            throw std::runtime_error ("npy: header without " + key);
        return header.substr (header.find (':', k) + 1);
    };
    {
        std::string d = field ("descr");
        const auto q0 = d.find ('\'');
        const auto q1 = d.find ('\'', q0 + 1);
        a.dtype = d.substr (q0 + 1, q1 - q0 - 1);
    }
    if (field ("fortran_order").find ("True") < field ("fortran_order").find (','))
        throw std::runtime_error ("npy: Fortran order is not supported");
    {
        std::string s = field ("shape");
        const auto p0 = s.find ('(');
        const auto p1 = s.find (')');
        std::string inside = s.substr (p0 + 1, p1 - p0 - 1);
        size_t pos = 0;
        while (pos < inside.size())
        {
            const auto comma = inside.find (',', pos);
            const std::string tok = inside.substr (pos, comma == std::string::npos ? std::string::npos : comma - pos);
            if (tok.find_first_not_of (" ") != std::string::npos)
                a.shape.push_back (static_cast<size_t> (std::stoull (tok)));
            if (comma == std::string::npos)
                break;
            pos = comma + 1;
        }
    }
    const int width = a.dtype.size() > 2 ? std::stoi (a.dtype.substr (2)) : 1;
    if (a.dtype[0] == '>')
        throw std::runtime_error ("npy: big-endian data are not supported");
    a.bytes.resize (a.count() * static_cast<size_t> (width));
    f.read (a.bytes.data(), static_cast<std::streamsize> (a.bytes.size()));
    if (! f)
        throw std::runtime_error ("npy: truncated data in " + path);
    return a;
}

// Writes a C-order .npy (format 1.0); T is float, double or uint8_t.
template <typename T>
void saveNpy (const std::string& path, const T* data, const std::vector<size_t>& shape)
{
    static_assert (sizeof (T) == 8 || sizeof (T) == 4 || sizeof (T) == 1, "npy: unsupported type");
    const char* descr = sizeof (T) == 8 ? "<f8" : sizeof (T) == 4 ? "<f4" : "|u1";
    std::string s = "(";
    size_t n = 1;
    for (auto d : shape)
    {
        s += std::to_string (d) + ",";
        n *= d;
    }
    if (shape.size() > 1)
        s.pop_back();
    s += ")";
    std::string header = std::string ("{'descr': '") + descr + "', 'fortran_order': False, 'shape': " + s + ", }";
    const size_t total = 10 + header.size() + 1;
    header.append ((64 - total % 64) % 64, ' ');
    header += '\n';
    std::ofstream f (path, std::ios::binary);
    if (! f)
        throw std::runtime_error ("npy: cannot write " + path);
    f.write ("\x93NUMPY\x01\x00", 8);
    const uint16_t len = static_cast<uint16_t> (header.size());
    f.write (reinterpret_cast<const char*> (&len), 2);
    f.write (header.data(), static_cast<std::streamsize> (header.size()));
    f.write (reinterpret_cast<const char*> (data), static_cast<std::streamsize> (n * sizeof (T)));
}
} // namespace ProbeViewer::Localization::Realtime
