#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <string>
#include <limits>
#include <random>

using namespace std;

constexpr double EPS = 1e-9;
constexpr double EPS2 = EPS * EPS;
constexpr double INF = 1e20;
const double PI = acos(-1.0);

inline double absd(double x) { return x < 0 ? -x : x; }

struct Pt
{
    double x, y;
    Pt() : x(0), y(0) {}
    Pt(double _x, double _y) : x(_x), y(_y) {}

    Pt operator+(const Pt &o) const { return Pt(x + o.x, y + o.y); }
    Pt operator-(const Pt &o) const { return Pt(x - o.x, y - o.y); }
    Pt operator-() const { return Pt(-x, -y); }
    Pt operator*(double k) const { return Pt(x * k, y * k); }
    Pt operator/(double k) const { return Pt(x / k, y / k); }
};

using Vec = Pt;

inline double cross(const Vec &a, const Vec &b) { return a.x * b.y - a.y * b.x; }
inline double dot(const Vec &a, const Vec &b) { return a.x * b.x + a.y * b.y; }
inline double len2(const Vec &v) { return dot(v, v); }

inline double dist2(const Pt &a, const Pt &b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

inline double signedArea(const vector<Pt> &p)
{
    double s = 0.0;
    int n = (int)p.size();
    for (int i = 0; i < n; ++i)
        s += cross(p[i], p[(i + 1) % n]);
    return s * 0.5;
}

inline bool pointOnSeg(const Pt &p, const Pt &a, const Pt &b)
{
    if (absd(cross(b - a, p - a)) > EPS)
        return false;
    return p.x >= min(a.x, b.x) - EPS && p.x <= max(a.x, b.x) + EPS &&
           p.y >= min(a.y, b.y) - EPS && p.y <= max(a.y, b.y) + EPS;
}

inline bool pointInPoly(const Pt &p, const vector<Pt> &poly)
{
    int n = (int)poly.size();
    if (n < 3)
        return false;

    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        const Pt &a = poly[i];
        const Pt &b = poly[j];

        if (pointOnSeg(p, a, b))
            return true;

        if ((a.y > p.y) != (b.y > p.y))
        {
            double x = (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x;
            if (p.x < x - EPS)
                inside = !inside;
        }
    }
    return inside;
}

struct BBox
{
    double minX, maxX, minY, maxY;
};

inline BBox getBBox(const vector<Pt> &poly)
{
    BBox b{INF, -INF, INF, -INF};
    for (const auto &p : poly)
    {
        b.minX = min(b.minX, p.x);
        b.maxX = max(b.maxX, p.x);
        b.minY = min(b.minY, p.y);
        b.maxY = max(b.maxY, p.y);
    }
    return b;
}

inline bool bboxOverlapFast(const BBox &b1, const BBox &b2, const Pt &t)
{
    return !(b1.maxX < b2.minX + t.x - EPS || b2.maxX + t.x < b1.minX - EPS ||
             b1.maxY < b2.minY + t.y - EPS || b2.maxY + t.y < b1.minY - EPS);
}

struct Polygon
{
    vector<Pt> v;
    int n;

    Polygon() : n(0) {}
    Polygon(const vector<Pt> &a) : v(a), n((int)a.size())
    {
        if (n > 0 && signedArea(v) < 0)
            reverse(v.begin(), v.end());
    }
    Polygon(vector<Pt> &&a) : v(move(a)), n((int)v.size())
    {
        if (n > 0 && signedArea(v) < 0)
            reverse(v.begin(), v.end());
    }

    double area() const { return signedArea(v); }

    vector<bool> convexMask() const
    {
        vector<bool> ok(n, true);
        if (n < 3)
            return ok;

        int s = 0;
        for (int i = 1; i < n; ++i)
        {
            if (v[i].y < v[s].y - EPS || (absd(v[i].y - v[s].y) < EPS && v[i].x < v[s].x))
                s = i;
        }

        Vec e1 = v[s] - v[(s - 1 + n) % n];
        Vec e2 = v[(s + 1) % n] - v[s];
        bool positive = cross(e1, e2) > EPS;

        for (int i = 0; i < n; ++i)
        {
            Vec a = v[i] - v[(i - 1 + n) % n];
            Vec b = v[(i + 1) % n] - v[i];
            double c = cross(a, b);
            if (absd(c) >= EPS)
                ok[i] = ((c > 0) == positive);
        }
        return ok;
    }

    struct BadRun
    {
        int l, r;
        Pt pl, pr;
        vector<Pt> seg;
    };

    vector<BadRun> findBadRuns() const
    {
        vector<BadRun> out;
        if (n < 3)
            return out;

        vector<bool> ok = convexMask();
        vector<int> bad;
        bad.reserve(n / 4);

        for (int i = 0; i < n; ++i)
            if (!ok[i])
                bad.push_back(i);

        if (bad.empty())
            return out;

        vector<pair<int, int>> runs;
        int L = bad[0], R = bad[0];
        for (size_t i = 1; i < bad.size(); ++i)
        {
            if (bad[i] == (R + 1) % n)
                R = bad[i];
            else
            {
                runs.push_back({L, R});
                L = R = bad[i];
            }
        }
        runs.push_back({L, R});

        for (auto [l, r] : runs)
        {
            int a = (l - 1 + n) % n;
            int b = (r + 1) % n;
            bool good = true;
            double sign = 0.0;

            for (int i = l;; i = (i + 1) % n)
            {
                double s = cross(v[b] - v[a], v[i] - v[a]);
                if (absd(s) > EPS)
                {
                    if (sign == 0.0)
                        sign = s;
                    else if ((sign > 0) != (s > 0))
                    {
                        good = false;
                        break;
                    }
                }
                if (i == r)
                    break;
            }

            if (!good)
            {
                a = 0;
                b = n - 1;
            }

            BadRun t;
            t.l = a;
            t.r = b;
            t.pl = v[a];
            t.pr = v[b];
            for (int i = a;; i = (i + 1) % n)
            {
                t.seg.push_back(v[i]);
                if (i == b)
                    break;
            }
            out.push_back(move(t));
        }

        return out;
    }

    Polygon convexHull() const
    {
        if (n <= 3)
            return *this;

        vector<Pt> a = v;
        sort(a.begin(), a.end(), [](const Pt &x, const Pt &y)
             {
                 if (x.x != y.x) return x.x < y.x;
                 return x.y < y.y; });

        a.erase(unique(a.begin(), a.end(), [](const Pt &x, const Pt &y)
                       { return dist2(x, y) < EPS2; }),
                a.end());

        if ((int)a.size() < 3)
            return *this;

        vector<Pt> h;
        h.reserve(a.size());

        for (const auto &p : a)
        {
            while (h.size() >= 2 &&
                   cross(h.back() - h[h.size() - 2], p - h.back()) > -EPS)
            {
                h.pop_back();
            }
            h.push_back(p);
        }

        size_t cut = h.size();
        for (int i = (int)a.size() - 2; i >= 0; --i)
        {
            while (h.size() > cut &&
                   cross(h.back() - h[h.size() - 2], a[i] - h.back()) > -EPS)
            {
                h.pop_back();
            }
            h.push_back(a[i]);
        }

        if (h.size() > 1 && dist2(h.front(), h.back()) < EPS2)
            h.pop_back();

        if (h.size() < 3)
            return *this;

        return Polygon(move(h));
    }

    Polygon repair() const
    {
        auto runs = findBadRuns();
        if (runs.empty())
            return *this;

        vector<char> keep(n, 1);
        for (const auto &r : runs)
        {
            for (int i = (r.l + 1) % n; i != r.r; i = (i + 1) % n)
                keep[i] = 0;
        }

        vector<Pt> r;
        r.reserve(n / 2);
        for (int i = 0; i < n; ++i)
            if (keep[i])
                r.push_back(v[i]);

        if (r.size() < 3)
            return convexHull();

        Polygon tmp(move(r));
        auto ok = tmp.convexMask();
        if (any_of(ok.begin(), ok.end(), [](bool x)
                   { return !x; }))
            return tmp.convexHull();

        return tmp;
    }
};

inline int lowestLeft(const Polygon &g)
{
    int id = 0;
    double my = g.v[0].y, mx = g.v[0].x;
    for (int i = 1; i < g.n; ++i)
    {
        double y = g.v[i].y, x = g.v[i].x;
        if (y < my - EPS || (absd(y - my) < EPS && x < mx - EPS))
        {
            id = i;
            my = y;
            mx = x;
        }
    }
    return id;
}

inline int highestRight(const Polygon &g)
{
    int id = 0;
    double my = g.v[0].y, mx = g.v[0].x;
    for (int i = 1; i < g.n; ++i)
    {
        double y = g.v[i].y, x = g.v[i].x;
        if (y > my + EPS || (absd(y - my) < EPS && x > mx + EPS))
        {
            id = i;
            my = y;
            mx = x;
        }
    }
    return id;
}

Polygon nfpConvex(const Polygon &a, const Polygon &b)
{
    if (a.n < 3 || b.n < 3)
        return Polygon();

    int sa = lowestLeft(a);
    int sb = highestRight(b);

    struct Edge
    {
        Vec v;
        double ang;
    };

    vector<Edge> ea, eb;
    ea.reserve(a.n);
    eb.reserve(b.n);

    for (int i = 0; i < a.n; ++i)
    {
        Vec e = a.v[(sa + i + 1) % a.n] - a.v[(sa + i) % a.n];
        if (len2(e) > EPS2)
            ea.push_back({e, atan2(e.y, e.x)});
    }

    for (int i = 0; i < b.n; ++i)
    {
        Vec e = b.v[(sb + i + 1) % b.n] - b.v[(sb + i) % b.n];
        e = -e;
        if (len2(e) > EPS2)
            eb.push_back({e, atan2(e.y, e.x)});
    }

    auto normAngle = [](double ang) -> double
    {
        while (ang < 0)
            ang += 2.0 * PI;
        while (ang >= 2.0 * PI)
            ang -= 2.0 * PI;
        return ang;
    };

    for (auto &x : ea)
        x.ang = normAngle(x.ang);
    for (auto &x : eb)
        x.ang = normAngle(x.ang);

    sort(ea.begin(), ea.end(), [](const Edge &l, const Edge &r)
         { return l.ang < r.ang; });
    sort(eb.begin(), eb.end(), [](const Edge &l, const Edge &r)
         { return l.ang < r.ang; });

    vector<Vec> mv;
    mv.reserve(ea.size() + eb.size());

    size_t i = 0, j = 0;
    while (i < ea.size() && j < eb.size())
    {
        if (ea[i].ang < eb[j].ang - EPS)
            mv.push_back(ea[i++].v);
        else if (eb[j].ang < ea[i].ang - EPS)
            mv.push_back(eb[j++].v);
        else
        {
            Vec s = ea[i].v + eb[j].v;
            if (len2(s) > EPS2)
                mv.push_back(s);
            ++i;
            ++j;
        }
    }
    while (i < ea.size())
        mv.push_back(ea[i++].v);
    while (j < eb.size())
        mv.push_back(eb[j++].v);

    if (mv.empty())
        return Polygon();

    Pt cur = a.v[sa] - b.v[sb];
    vector<Pt> h;
    h.reserve(mv.size() + 1);
    h.push_back(cur);

    for (const auto &e : mv)
    {
        cur = cur + e;
        h.push_back(cur);
    }

    if (h.size() > 2 && dist2(h.front(), h.back()) < EPS2)
        h.pop_back();

    if (h.size() < 3)
        return Polygon();

    double s = 0.0;
    for (size_t k = 0; k < h.size(); ++k)
        s += cross(h[k], h[(k + 1) % h.size()]);

    if (s < 0)
        reverse(h.begin(), h.end());

    return Polygon(move(h));
}

Polygon nfpAll(const Polygon &a, const Polygon &b)
{
    auto ca = a.convexMask();
    auto cb = b.convexMask();

    bool fa = all_of(ca.begin(), ca.end(), [](bool x)
                     { return x; });
    bool fb = all_of(cb.begin(), cb.end(), [](bool x)
                     { return x; });

    if (fa && fb)
        return nfpConvex(a, b);

    Polygon ra = a.repair();
    Polygon rb = b.repair();
    Polygon r = nfpConvex(ra, rb);

    if (r.n < 3)
        return Polygon();

    return r;
}

int kDirectionCount = 32;
vector<Pt> g_escapeDirs; // 存的是“平移空间”的候选边界点
BBox bbox1, bbox2;
Pt g_ref; // polygon2 的参考点

vector<Pt> BuildDirections(const Polygon &polygon1, const Polygon &polygon2)
{
    vector<Pt> dirs;
    dirs.reserve(kDirectionCount);

    // 固定随机种子，保证每次运行结果一致
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 2.0 * PI);

    // 随机一个起始角度，再均匀旋转一圈
    double base = dist(rng);

    for (int i = 0; i < kDirectionCount; ++i)
    {
        double ang = base + 2.0 * PI * i / kDirectionCount;
        dirs.emplace_back(cos(ang), sin(ang));
    }

    return dirs;
}

// vector<Pt> BuildDirections(const Polygon &polygon1, const Polygon &polygon2)
// {
//     struct Cand
//     {
//         double score; // 用边长平方衡量这条法向的重要性
//         Pt dir;
//         Cand() : score(-1.0), dir(0, 0) {}
//     };

//     const int BIN_CNT = 128;
//     const double TWO_PI = 2.0 * PI;
//     vector<Cand> bins(BIN_CNT);

//     auto normAngle = [&](double ang) -> double
//     {
//         while (ang < 0)
//             ang += TWO_PI;
//         while (ang >= TWO_PI)
//             ang -= TWO_PI;
//         return ang;
//     };

//     auto addPoly = [&](const Polygon &poly)
//     {
//         int n = poly.n;
//         for (int i = 0; i < n; ++i)
//         {
//             Pt a = poly.v[i];
//             Pt b = poly.v[(i + 1) % n];
//             Pt e = b - a;

//             double e2 = len2(e);
//             if (e2 < 1e-18)
//                 continue;

//             Pt nrm(-e.y, e.x);
//             double L = sqrt(e2);
//             nrm = nrm * (1.0 / L);

//             double ang = atan2(nrm.y, nrm.x);
//             ang = normAngle(ang);

//             int bin = (int)(ang / TWO_PI * BIN_CNT);
//             if (bin < 0)
//                 bin = 0;
//             if (bin >= BIN_CNT)
//                 bin = BIN_CNT - 1;

//             if (e2 > bins[bin].score)
//             {
//                 bins[bin].score = e2;
//                 bins[bin].dir = nrm;
//             }
//         }
//     };

//     addPoly(polygon1);
//     addPoly(polygon2);

//     vector<Pt> cand;
//     cand.reserve(BIN_CNT);
//     for (int i = 0; i < BIN_CNT; ++i)
//         if (bins[i].score > 0)
//             cand.push_back(bins[i].dir);

//     vector<Pt> dirs;
//     dirs.reserve(kDirectionCount);

//     if ((int)cand.size() >= kDirectionCount)
//     {
//         for (int i = 0; i < kDirectionCount; ++i)
//         {
//             int idx = (long long)i * (int)cand.size() / kDirectionCount;
//             if (idx >= (int)cand.size())
//                 idx = (int)cand.size() - 1;
//             dirs.push_back(cand[idx]);
//         }
//     }
//     else
//     {
//         for (const auto &v : cand)
//             dirs.push_back(v);

//         for (int i = 0; (int)dirs.size() < kDirectionCount; ++i)
//         {
//             double ang = TWO_PI * i / kDirectionCount;
//             dirs.emplace_back(cos(ang), sin(ang));
//         }
//     }

//     return dirs;
// }

// 从 NFP 上沿 dir 射线求第一个边界交点，返回“绝对坐标”的边界点
Pt GenSolution(const Polygon &nfp, const Pt &dir)
{
    if (nfp.n < 3)
        return Pt(0, 0);

    Pt d = dir;
    double L = sqrt(len2(d));
    if (L < EPS)
        return Pt(0, 0);
    d = d * (1.0 / L);

    // 取一个内部点：凸多边形顶点平均值一定在内部
    Pt center(0, 0);
    for (const auto &p : nfp.v)
        center = center + p;
    center = center / (double)nfp.n;

    double bestT = INF;
    Pt bestP = center;

    for (int i = 0; i < nfp.n; ++i)
    {
        Pt a = nfp.v[i];
        Pt b = nfp.v[(i + 1) % nfp.n];
        Pt s = b - a;

        double den = cross(d, s);
        if (absd(den) < EPS)
            continue;

        // center + t*d = a + u*s
        double t = cross(a - center, s) / den;
        double u = cross(a - center, d) / den;

        if (t >= -EPS && u >= -EPS && u <= 1.0 + EPS)
        {
            if (t < bestT)
            {
                bestT = t;
                bestP = center + d * t;
            }
        }
    }

    return bestP;
}

void PreProcess(const Polygon &polygon1, const Polygon &polygon2, const Polygon &nfp)
{
    bbox1 = getBBox(polygon1.v);
    bbox2 = getBBox(polygon2.v);

    int base = polygon1.n + polygon2.n;
    kDirectionCount = 25;

    g_ref = polygon2.v.empty() ? Pt(0, 0) : polygon2.v[0];

    vector<Pt> dirs = BuildDirections(polygon1, polygon2);
    g_escapeDirs.reserve(kDirectionCount);

    for (int i = 0; i < kDirectionCount; ++i)
    {
        Pt dir = dirs[i];
        Pt q = GenSolution(nfp, dir);      // 绝对坐标上的边界点
        g_escapeDirs.push_back(q - g_ref); // 转成平移空间
    }
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    int n1, n2;
    cin >> n1 >> n2;

    vector<Pt> aPts(n1), bPts(n2);
    for (int i = 0; i < n1; ++i)
        cin >> aPts[i].x >> aPts[i].y;
    for (int i = 0; i < n2; ++i)
        cin >> bPts[i].x >> bPts[i].y;

    string ok;
    cin >> ok;

    Polygon polygon1(move(aPts));
    Polygon polygon2(move(bPts));

    Polygon nfp = nfpAll(polygon1, polygon2);

    PreProcess(polygon1, polygon2, nfp);

    cout << "OK\n";
    cout.flush();

    int m;
    cin >> m;

    vector<Pt> testCases(m);
    for (int i = 0; i < m; ++i)
        cin >> testCases[i].x >> testCases[i].y;

    cin >> ok;

    cout << m << '\n';

    string buffer;
    buffer.reserve(m * 32); // 预留空间，避免频繁扩容

    for (int i = 0; i < m; ++i)
    {
        const Pt &t = testCases[i];

        if (!bboxOverlapFast(bbox1, bbox2, t))
        {
            buffer += "0.00000 0.00000\n";
            continue;
        }

        Pt res;
        double best = numeric_limits<double>::infinity();

        for (int j = 0; j < kDirectionCount; ++j)
        {
            Pt tmp = g_escapeDirs[j] - t;
            double d = len2(tmp);
            if (d < best)
            {
                best = d;
                res = tmp;
            }
        }

        char line[64];
        snprintf(line, sizeof(line), "%.5f %.5f\n", res.x, res.y);
        buffer += line;
    }

    // 一次性输出
    cout << buffer;

    cout << "OK\n";
    cout.flush();
    return 0;
}