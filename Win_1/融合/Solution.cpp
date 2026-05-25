#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include <thread>
#include <chrono>

struct Vector2D;
struct Polygon;

using namespace std;

const double EPS = 1e-6;
const double eps = 1e-9;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double EPS_SQ = EPS * EPS;
constexpr double INF_D = 1e20;

struct Vector2D
{
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}
    Vector2D operator-(const Vector2D &other) const { return Vector2D(x - other.x, y - other.y); }
    Vector2D operator+(const Vector2D &other) const { return Vector2D(x + other.x, y + other.y); }
    Vector2D operator*(double scalar) const { return Vector2D(x * scalar, y * scalar); }
    Vector2D operator/(double scalar) const { return Vector2D(x / scalar, y / scalar); }
    Vector2D &operator+=(const Vector2D &other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vector2D &operator-=(const Vector2D &other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    double Dot(const Vector2D &other) const { return x * other.x + y * other.y; }
    double Length() const { return std::sqrt(x * x + y * y); }
    Vector2D Normalize() const
    {
        double len = Length();
        if (len == 0)
            return *this;
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

struct Polygon
{
    std::vector<Vector2D> vertices;

    Polygon() = default;
    Polygon(std::initializer_list<Vector2D> vts) : vertices(vts) {}
    Polygon(std::vector<Vector2D> vts) : vertices(std::move(vts)) {}

    Vector2D GetCenter() const
    {
        Vector2D center(0, 0);
        if (vertices.empty())
            return center;
        for (const auto &v : vertices)
            center = center + v;
        return center * (1.0 / vertices.size());
    }
    void MoveByVec(const Vector2D &vec)
    {
        for (auto &v : vertices)
            v = v + vec;
    }
};

// ------------------------
// 基础几何
// ------------------------

inline double fast_abs(double x) { return x < 0 ? -x : x; }

double Cross(const Vector2D &a, const Vector2D &b) { return a.x * b.y - a.y * b.x; }
double Cross(const Vector2D &o, const Vector2D &a, const Vector2D &b) { return Cross(a - o, b - o); }
inline double Dot(const Vector2D &a, const Vector2D &b) { return a.x * b.x + a.y * b.y; }
inline double Length2(const Vector2D &v) { return v.x * v.x + v.y * v.y; }
inline double Distance2(const Vector2D &a, const Vector2D &b) { return Length2(a - b); }

vector<Vector2D> ConvexHull(vector<Vector2D> pts)
{
    int n = (int)pts.size();
    if (n <= 2)
        return pts;
    sort(pts.begin(), pts.end(), [](const Vector2D &a, const Vector2D &b)
         {
        if (fabs(a.x - b.x) > EPS) return a.x < b.x;
        return a.y < b.y; });
    pts.erase(unique(pts.begin(), pts.end(), [](const Vector2D &a, const Vector2D &b)
                     { return Distance2(a, b) < EPS_SQ; }),
              pts.end());
    if ((int)pts.size() <= 2)
        return pts;

    vector<Vector2D> hull;
    for (int i = 0; i < (int)pts.size(); ++i)
    {
        while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= EPS)
            hull.pop_back();
        hull.push_back(pts[i]);
    }
    int t = (int)hull.size();
    for (int i = (int)pts.size() - 2; i >= 0; --i)
    {
        while ((int)hull.size() > t && Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= EPS)
            hull.pop_back();
        hull.push_back(pts[i]);
    }
    if (!hull.empty())
        hull.pop_back();
    return hull;
}

inline double polygonArea(const vector<Vector2D> &poly)
{
    double res = 0;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i)
        res += Cross(poly[i], poly[(i + 1) % n]);
    return res / 2.0;
}

inline double polygonArea(const Polygon &poly)
{
    return polygonArea(poly.vertices);
}

inline bool onSegment(const Vector2D &A, const Vector2D &B, const Vector2D &P)
{
    if (std::fabs(Cross(A, B, P)) > eps)
        return false;
    double minX = std::min(A.x, B.x) - eps;
    double maxX = std::max(A.x, B.x) + eps;
    double minY = std::min(A.y, B.y) - eps;
    double maxY = std::max(A.y, B.y) + eps;
    return (P.x >= minX && P.x <= maxX && P.y >= minY && P.y <= maxY);
}

bool PointInPolygonStrict(const Vector2D &p, const Polygon &poly)
{
    const auto &v = poly.vertices;
    int n = (int)v.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        const Vector2D &a = v[i];
        const Vector2D &b = v[j];
        if (onSegment(a, b, p))
            return false;
        if ((a.y > p.y) != (b.y > p.y))
        {
            double x = a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y);
            if (x > p.x)
                inside = !inside;
        }
    }
    return inside;
}

int SegmentsStrictIntersect(const Vector2D &a, const Vector2D &b, const Vector2D &c, const Vector2D &d)
{
    double d1 = Cross(a, b, c);
    double d2 = Cross(a, b, d);
    double d3 = Cross(c, d, a);
    double d4 = Cross(c, d, b);
    if (d1 * d2 < -EPS && d3 * d4 < -EPS)
        return 5;
    if (onSegment(c, d, a))
        return 1;
    if (onSegment(c, d, b))
        return 2;
    if (onSegment(a, b, c))
        return 3;
    if (onSegment(a, b, d))
        return 4;
    return 0;
}

struct Projection
{
    double min, max;
};

Projection ProjectPolygon(const Polygon &poly, const Vector2D &axis)
{
    double minProj = poly.vertices[0].Dot(axis);
    double maxProj = minProj;
    for (size_t i = 1; i < poly.vertices.size(); ++i)
    {
        double proj = poly.vertices[i].Dot(axis);
        if (proj < minProj)
            minProj = proj;
        if (proj > maxProj)
            maxProj = proj;
    }
    return {minProj, maxProj};
}

double calOverlapWithDir(const Projection &projA, const Projection &projB)
{
    if (projA.max <= projB.min || projB.max <= projA.min)
        return 0.0;
    double d1 = projA.max - projB.min;
    double d2 = projB.max - projA.min;
    return (d1 < d2) ? d1 : d2;
}

bool bboxOverlap(const Polygon &P, const Polygon &Q)
{
    if (P.vertices.empty() || Q.vertices.empty())
        return false;
    double minXP = P.vertices[0].x, maxXP = P.vertices[0].x, minYP = P.vertices[0].y, maxYP = P.vertices[0].y;
    double minXQ = Q.vertices[0].x, maxXQ = Q.vertices[0].x, minYQ = Q.vertices[0].y, maxYQ = Q.vertices[0].y;
    for (const auto &v : P.vertices)
    {
        minXP = std::min(minXP, v.x);
        maxXP = std::max(maxXP, v.x);
        minYP = std::min(minYP, v.y);
        maxYP = std::max(maxYP, v.y);
    }
    for (const auto &v : Q.vertices)
    {
        minXQ = std::min(minXQ, v.x);
        maxXQ = std::max(maxXQ, v.x);
        minYQ = std::min(minYQ, v.y);
        maxYQ = std::max(maxYQ, v.y);
    }
    return !(maxXP + eps < minXQ || maxXQ + eps < minXP || maxYP + eps < minYQ || maxYQ + eps < minYP);
}

Vector2D movepoint(Vector2D node, Vector2D vec) { return node + vec; }

bool check(const Polygon &P, const Polygon &Q, const Vector2D &vec)
{
    Polygon moved = Q;
    moved.MoveByVec(vec);

    if (!bboxOverlap(P, moved))
        return true;

    const auto &vP = P.vertices;
    const auto &vQ = moved.vertices;
    int nP = (int)vP.size(), nQ = (int)vQ.size();
    for (int i = 0; i < nP; ++i)
    {
        const Vector2D &a = vP[i];
        const Vector2D &b = vP[(i + 1) % nP];
        for (int j = 0; j < nQ; ++j)
        {
            const Vector2D &c = vQ[j];
            const Vector2D &d = vQ[(j + 1) % nQ];
            int op = SegmentsStrictIntersect(a, b, c, d);
            if (op)
                return false;
        }
    }
    for (const auto &q : vQ)
        if (PointInPolygonStrict(q, P))
            return false;
    for (const auto &p : vP)
        if (PointInPolygonStrict(p, moved))
            return false;
    return true;
}

// ------------------------
// NFP 相关模块
// ------------------------

bool isConvexPolygon(const Polygon &poly)
{
    int n = (int)poly.vertices.size();
    if (n < 3)
        return false;
    int sign = 0;
    for (int i = 0; i < n; ++i)
    {
        const Vector2D &a = poly.vertices[i];
        const Vector2D &b = poly.vertices[(i + 1) % n];
        const Vector2D &c = poly.vertices[(i + 2) % n];
        double cr = Cross(a, b, c);
        if (fabs(cr) < EPS)
            continue;
        int cur = (cr > 0) ? 1 : -1;
        if (sign == 0)
            sign = cur;
        else if (sign != cur)
            return false;
    }
    return true;
}

Polygon normalizeCCW(Polygon poly)
{
    if (poly.vertices.size() >= 3 && polygonArea(poly) < 0)
        std::reverse(poly.vertices.begin(), poly.vertices.end());
    return poly;
}

int findLowestLeftmost(const Polygon &poly)
{
    int idx = 0;
    double minY = poly.vertices[0].y;
    double minX = poly.vertices[0].x;
    for (int i = 1; i < (int)poly.vertices.size(); ++i)
    {
        double y = poly.vertices[i].y;
        double x = poly.vertices[i].x;
        if (y < minY - EPS || (fast_abs(y - minY) < EPS && x < minX - EPS))
        {
            idx = i;
            minY = y;
            minX = x;
        }
    }
    return idx;
}

int findHighestRightmost(const Polygon &poly)
{
    int idx = 0;
    double maxY = poly.vertices[0].y;
    double maxX = poly.vertices[0].x;
    for (int i = 1; i < (int)poly.vertices.size(); ++i)
    {
        double y = poly.vertices[i].y;
        double x = poly.vertices[i].x;
        if (y > maxY + EPS || (fast_abs(y - maxY) < EPS && x > maxX + EPS))
        {
            idx = i;
            maxY = y;
            maxX = x;
        }
    }
    return idx;
}

Polygon nfpMoving(const Polygon &A, const Polygon &B)
{
    if (A.vertices.size() < 3 || B.vertices.size() < 3)
        return Polygon();

    Polygon a = normalizeCCW(A);
    Polygon b = normalizeCCW(B);
    int startA = findLowestLeftmost(a);
    int startB = findHighestRightmost(b);

    vector<std::pair<Vector2D, double>> edgesA, edgesB;
    edgesA.reserve(a.vertices.size());
    edgesB.reserve(b.vertices.size());

    for (int i = 0; i < (int)a.vertices.size(); ++i)
    {
        Vector2D e = a.vertices[(startA + i + 1) % a.vertices.size()] - a.vertices[(startA + i) % a.vertices.size()];
        if (Length2(e) > EPS_SQ)
            edgesA.push_back({e, std::atan2(e.y, e.x)});
    }
    for (int i = 0; i < (int)b.vertices.size(); ++i)
    {
        Vector2D e = b.vertices[(startB + i + 1) % b.vertices.size()] - b.vertices[(startB + i) % b.vertices.size()];
        e = e * -1.0;
        if (Length2(e) > EPS_SQ)
            edgesB.push_back({e, std::atan2(e.y, e.x)});
    }

    sort(edgesA.begin(), edgesA.end(), [](const auto &lhs, const auto &rhs)
         { return lhs.second < rhs.second; });
    sort(edgesB.begin(), edgesB.end(), [](const auto &lhs, const auto &rhs)
         { return lhs.second < rhs.second; });

    vector<Vector2D> merged;
    merged.reserve(edgesA.size() + edgesB.size());
    size_t i = 0, j = 0;
    while (i < edgesA.size() && j < edgesB.size())
    {
        if (edgesA[i].second < edgesB[j].second - EPS)
            merged.push_back(edgesA[i++].first);
        else if (edgesB[j].second < edgesA[i].second - EPS)
            merged.push_back(edgesB[j++].first);
        else
        {
            Vector2D m = edgesA[i].first + edgesB[j].first;
            if (Length2(m) > EPS_SQ)
                merged.push_back(m);
            ++i;
            ++j;
        }
    }
    while (i < edgesA.size())
        merged.push_back(edgesA[i++].first);
    while (j < edgesB.size())
        merged.push_back(edgesB[j++].first);

    if (merged.empty())
        return Polygon();

    Vector2D cur = a.vertices[startA] - b.vertices[startB];
    vector<Vector2D> nfp;
    nfp.reserve(merged.size() + 1);
    nfp.push_back(cur);
    for (const auto &e : merged)
    {
        cur = cur + e;
        nfp.push_back(cur);
    }

    if (nfp.size() > 2 && Distance2(nfp.front(), nfp.back()) < EPS_SQ)
        nfp.pop_back();
    if (nfp.size() < 3)
        return Polygon();
    if (polygonArea(nfp) < 0)
        std::reverse(nfp.begin(), nfp.end());
    return Polygon(std::move(nfp));
}

struct ConcaveSet
{
    int startIdx, endIdx;
    Vector2D replaceStart, replaceEnd;
    vector<Vector2D> localContour;
};

vector<bool> getConvexityMask(const Polygon &poly)
{
    int n = (int)poly.vertices.size();
    vector<bool> conv(n, true);
    if (n < 3)
        return conv;

    int minYIdx = 0;
    double minY = poly.vertices[0].y;
    for (int i = 1; i < n; ++i)
    {
        if (poly.vertices[i].y < minY - EPS || (fast_abs(poly.vertices[i].y - minY) < EPS && poly.vertices[i].x < poly.vertices[minYIdx].x))
        {
            minYIdx = i;
            minY = poly.vertices[i].y;
        }
    }

    Vector2D ePrev = poly.vertices[minYIdx] - poly.vertices[(minYIdx - 1 + n) % n];
    Vector2D eNext = poly.vertices[(minYIdx + 1) % n] - poly.vertices[minYIdx];
    bool expectPositive = Cross(ePrev, eNext) > EPS;

    for (int i = 0; i < n; ++i)
    {
        Vector2D e1 = poly.vertices[i] - poly.vertices[(i - 1 + n) % n];
        Vector2D e2 = poly.vertices[(i + 1) % n] - poly.vertices[i];
        double cr = Cross(e1, e2);
        if (fast_abs(cr) >= EPS)
            conv[i] = (cr > 0) == expectPositive;
    }
    return conv;
}

vector<ConcaveSet> extractConcaveSets(const Polygon &poly)
{
    vector<ConcaveSet> sets;
    int n = (int)poly.vertices.size();
    if (n < 3)
        return sets;

    vector<bool> conv = getConvexityMask(poly);
    vector<int> bad;
    for (int i = 0; i < n; ++i)
        if (!conv[i])
            bad.push_back(i);
    if (bad.empty())
        return sets;

    vector<std::pair<int, int>> groups;
    int gs = bad[0], ge = bad[0];
    for (size_t i = 1; i < bad.size(); ++i)
    {
        if (bad[i] == (ge + 1) % n)
            ge = bad[i];
        else
        {
            groups.push_back({gs, ge});
            gs = ge = bad[i];
        }
    }
    groups.push_back({gs, ge});

    for (const auto &g : groups)
    {
        int a = (g.first - 1 + n) % n, b = (g.second + 1) % n;
        bool valid = true;
        double expSide = 0;
        for (int i = g.first;; i = (i + 1) % n)
        {
            double side = Cross(poly.vertices[a], poly.vertices[b], poly.vertices[i]);
            if (fast_abs(side) > EPS)
            {
                if (expSide == 0)
                    expSide = side;
                else if ((expSide > 0) != (side > 0))
                {
                    valid = false;
                    break;
                }
            }
            if (i == g.second)
                break;
        }
        if (!valid)
        {
            a = 0;
            b = n - 1;
        }

        ConcaveSet cs;
        cs.startIdx = a;
        cs.endIdx = b;
        cs.replaceStart = poly.vertices[a];
        cs.replaceEnd = poly.vertices[b];
        for (int i = a;; i = (i + 1) % n)
        {
            cs.localContour.push_back(poly.vertices[i]);
            if (i == b)
                break;
        }
        sets.push_back(std::move(cs));
    }
    return sets;
}

Polygon convexifyPolygon(const Polygon &poly)
{
    auto sets = extractConcaveSets(poly);
    if (sets.empty())
        return poly;

    int n = (int)poly.vertices.size();
    vector<bool> keep(n, true);
    for (const auto &cs : sets)
    {
        for (int i = (cs.startIdx + 1) % n; i != cs.endIdx; i = (i + 1) % n)
            keep[i] = false;
    }

    vector<Vector2D> nv;
    for (int i = 0; i < n; ++i)
        if (keep[i])
            nv.push_back(poly.vertices[i]);
    if (nv.size() < 3)
    {
        Polygon hull;
        hull.vertices = ConvexHull(poly.vertices);
        return normalizeCCW(hull);
    }

    Polygon res;
    res.vertices = std::move(nv);
    if (!isConvexPolygon(res))
    {
        Polygon hull;
        hull.vertices = ConvexHull(res.vertices);
        return normalizeCCW(hull);
    }
    return normalizeCCW(res);
}

Polygon toConvexHullPolygon(const Polygon &poly)
{
    Polygon hull;
    hull.vertices = ConvexHull(poly.vertices);
    return normalizeCCW(hull);
}

Polygon minkowskiDiff(const Polygon &A, const Polygon &B)
{
    vector<Vector2D> pts;
    pts.reserve(A.vertices.size() * B.vertices.size());
    for (int i = 0; i < (int)A.vertices.size(); ++i)
        for (int j = 0; j < (int)B.vertices.size(); ++j)
            pts.emplace_back(A.vertices[i].x - B.vertices[j].x, A.vertices[i].y - B.vertices[j].y);

    if (pts.size() < 3)
        return Polygon();
    sort(pts.begin(), pts.end(), [](const Vector2D &a, const Vector2D &b)
         {
        if (fabs(a.x - b.x) > EPS) return a.x < b.x;
        return a.y < b.y; });
    pts.erase(unique(pts.begin(), pts.end(), [](const Vector2D &a, const Vector2D &b)
                     { return Distance2(a, b) < EPS_SQ; }),
              pts.end());
    if (pts.size() < 3)
        return Polygon();

    Vector2D minPt = pts[0];
    for (const auto &p : pts)
        if (p.y < minPt.y - EPS || (fast_abs(p.y - minPt.y) < EPS && p.x < minPt.x - EPS))
            minPt = p;

    sort(pts.begin(), pts.end(), [&minPt](const Vector2D &a, const Vector2D &b)
         { return std::atan2(a.y - minPt.y, a.x - minPt.x) < std::atan2(b.y - minPt.y, b.x - minPt.x); });

    vector<Vector2D> hull;
    hull.reserve(pts.size());
    for (const auto &p : pts)
    {
        while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), p) < EPS)
            hull.pop_back();
        hull.push_back(p);
    }
    if (hull.size() > 1 && Distance2(hull.back(), hull[0]) < EPS_SQ)
        hull.pop_back();
    if (hull.size() < 3)
        return Polygon();
    if (polygonArea(hull) < 0)
        std::reverse(hull.begin(), hull.end());
    return Polygon(std::move(hull));
}

struct TrajSeg
{
    Vector2D start, end;
    double angle;
    double dx, dy;
    TrajSeg() : angle(0), dx(0), dy(0) {}
    TrajSeg(const Vector2D &s, const Vector2D &e) : start(s), end(e), dx(e.x - s.x), dy(e.y - s.y)
    {
        angle = std::atan2(dy, dx);
    }
    double len2() const { return dx * dx + dy * dy; }
};

inline bool lineIntWithT(const Vector2D &p1, const Vector2D &r, const Vector2D &q1, const Vector2D &s,
                         double &t_out, double &u_out, Vector2D &res)
{
    double rxs = Cross(r, s);
    if (fast_abs(rxs) < EPS)
        return false;
    Vector2D qp = q1 - p1;
    t_out = Cross(qp, s) / rxs;
    u_out = Cross(qp, r) / rxs;
    res = p1 + r * t_out;
    return true;
}

inline bool hasInt(const TrajSeg &t, const vector<TrajSeg> &all, int idx)
{
    for (size_t i = 0; i < all.size(); ++i)
    {
        if ((int)i == idx)
            continue;
        const auto &seg = all[i];
        Vector2D inter;
        double t1, t2;
        if (lineIntWithT(t.start, Vector2D(t.dx, t.dy), seg.start, Vector2D(seg.dx, seg.dy), t1, t2, inter))
        {
            if (t1 >= -EPS && t1 <= 1.0 + EPS && t2 >= -EPS && t2 <= 1.0 + EPS)
                return true;
        }
    }
    return false;
}

Polygon nfpTrajectoryFull(const Polygon &A, const Polygon &B)
{
    const int An = (int)A.vertices.size();
    const int Bn = (int)B.vertices.size();
    vector<TrajSeg> all;
    all.reserve(An * Bn * 2);
    if (An < 3 || Bn < 3)
        return Polygon();

    Vector2D refB = B.vertices[0];

    for (int i = 0; i < An; ++i)
    {
        Vector2D a1 = A.vertices[i], a2 = A.vertices[(i + 1) % An];
        double dx_a = a2.x - a1.x, dy_a = a2.y - a1.y;
        for (int j = 0; j < Bn; ++j)
        {
            Vector2D bj = B.vertices[j];
            Vector2D ts(a1.x - bj.x + refB.x, a1.y - bj.y + refB.y);
            Vector2D te(a2.x - bj.x + refB.x, a2.y - bj.y + refB.y);
            if (dx_a * dx_a + dy_a * dy_a > EPS_SQ)
                all.emplace_back(ts, te);
        }
    }
    for (int i = 0; i < Bn; ++i)
    {
        Vector2D b1 = B.vertices[i], b2 = B.vertices[(i + 1) % Bn];
        double dx_b = b2.x - b1.x, dy_b = b2.y - b1.y;
        for (int j = 0; j < An; ++j)
        {
            Vector2D aj = A.vertices[j];
            Vector2D ts(aj.x - b1.x + refB.x, aj.y - b1.y + refB.y);
            Vector2D te(aj.x - b2.x + refB.x, aj.y - b2.y + refB.y);
            if (dx_b * dx_b + dy_b * dy_b > EPS_SQ)
                all.emplace_back(ts, te);
        }
    }

    if (all.empty())
        return Polygon();

    int startIdx = -1;
    const size_t allSize = all.size();
    for (size_t i = 0; i < allSize; ++i)
    {
        if (hasInt(all[i], all, (int)i))
        {
            startIdx = (int)i;
            break;
        }
    }
    if (startIdx == -1)
        startIdx = 0;

    vector<bool> used(allSize, false);
    vector<Vector2D> nfp;
    int curIdx = startIdx;
    Vector2D curPos = all[startIdx].start;

    do
    {
        nfp.push_back(curPos);
        used[curIdx] = true;
        int nextIdx = -1;
        for (size_t i = 0; i < allSize; ++i)
        {
            if (used[i])
                continue;
            if (Distance2(all[i].start, curPos) < EPS_SQ)
            {
                nextIdx = (int)i;
                break;
            }
        }
        if (nextIdx == -1)
            break;
        curPos = all[nextIdx].end;
        curIdx = nextIdx;
    } while (curIdx != startIdx && nfp.size() < allSize * 2);

    if (nfp.size() < 3)
        return Polygon();
    if (polygonArea(nfp) < 0)
        std::reverse(nfp.begin(), nfp.end());
    return Polygon(std::move(nfp));
}

Polygon nfpCombined(const Polygon &A, const Polygon &B)
{
    bool isA = isConvexPolygon(A);
    bool isB = isConvexPolygon(B);
    if (isA && isB)
        return nfpMoving(A, B);

    Polygon Ac = convexifyPolygon(A);
    Polygon Bc = convexifyPolygon(B);
    Polygon nfpBase = nfpMoving(Ac, Bc);
    if (nfpBase.vertices.size() < 3)
        return Polygon();

    auto sA = extractConcaveSets(A);
    auto sB = extractConcaveSets(B);
    if (sA.empty() && sB.empty())
        return nfpBase;
    if ((int)sA.size() + (int)sB.size() <= 3)
        return nfpTrajectoryFull(A, B);
    return minkowskiDiff(A, B);
}

struct DistanceResult
{
    Vector2D nearest;
    double distance;
    bool inside;
};

DistanceResult pointToPolygon(const Vector2D &p, const Polygon &poly)
{
    DistanceResult res;
    res.inside = PointInPolygonStrict(p, poly);
    res.distance = INF_D;
    res.nearest = p;

    const int n = (int)poly.vertices.size();
    if (n < 3)
        return res;

    double minDist2 = INF_D;
    double bestDotY = -INF_D;
    double bestDotX = -INF_D;

    for (int i = 0; i < n; ++i)
    {
        Vector2D near;
        const Vector2D &a = poly.vertices[i];
        const Vector2D &b = poly.vertices[(i + 1) % n];
        Vector2D ab = b - a, ap = p - a;
        double ab2 = Dot(ab, ab);
        double d2;

        if (ab2 < EPS_SQ)
        {
            near = a;
            d2 = Dot(ap, ap);
        }
        else
        {
            double t = Dot(ap, ab) / ab2;
            if (t < 0)
            {
                near = a;
                d2 = Dot(ap, ap);
            }
            else if (t > 1)
            {
                near = b;
                Vector2D bp = p - b;
                d2 = Dot(bp, bp);
            }
            else
            {
                near = a + ab * t;
                Vector2D np = p - near;
                d2 = Dot(np, np);
            }
        }

        Vector2D vec = near - p;
        double dotY = vec.y;
        double dotX = vec.x;

        if (d2 < minDist2 - EPS_SQ)
        {
            minDist2 = d2;
            res.nearest = near;
            bestDotY = dotY;
            bestDotX = dotX;
        }
        else if (fast_abs(d2 - minDist2) < EPS_SQ)
        {
            if (dotY > bestDotY + EPS)
            {
                minDist2 = d2;
                res.nearest = near;
                bestDotY = dotY;
                bestDotX = dotX;
            }
            else if (fast_abs(dotY - bestDotY) < EPS && dotX > bestDotX + EPS)
            {
                minDist2 = d2;
                res.nearest = near;
                bestDotY = dotY;
                bestDotX = dotX;
            }
        }
    }

    res.distance = std::sqrt(minDist2);
    if (minDist2 < EPS_SQ)
        res.nearest = p;
    return res;
}

// ------------------------
// 原始方法：方向采样 + 二分
// ------------------------

int n1 = 0, n2 = 0, m = 0;
bool f1 = false, f2 = false;
Polygon polygon1;
Polygon polygon2;
Polygon polygon3, polygon4;
vector<Vector2D> testCases;
vector<Vector2D> axes, allaxes;

// 预计算的方向及对应的最小分离距离
std::vector<Vector2D> g_escapeDirs;
std::vector<double> g_escapeDist;

bool useNfpMode = false;
Polygon g_nfpPoly;
Vector2D g_refB;

std::vector<Vector2D> BuildUniformDirections(int cnt, double initAngle = 0.0)
{
    std::vector<Vector2D> out;
    out.reserve(cnt);
    const double step = 2.0 * PI / cnt;
    for (int i = 0; i < cnt; ++i)
    {
        double theta = initAngle + i * step;
        out.emplace_back(std::cos(theta), std::sin(theta));
    }
    return out;
}

Vector2D movepoint(Vector2D node, Vector2D vec)
{
    return node + vec;
}

double GenSolution(const Vector2D &dir)
{
    Vector2D unit = dir.Normalize();
    if (Length2(unit) < EPS_SQ)
        return 0.0;

    if (check(polygon1, polygon2, Vector2D(0, 0)))
        return 0.0;

    double lo = 0.0, hi = 1.0;
    while (hi < 1e6 && !check(polygon1, polygon2, unit * hi))
        hi *= 2.0;
    if (hi >= 1e6 && !check(polygon1, polygon2, unit * hi))
        return hi;

    for (int iter = 0; iter < 70; ++iter)
    {
        double mid = (lo + hi) * 0.5;
        if (check(polygon1, polygon2, unit * mid))
            hi = mid;
        else
            lo = mid;
    }
    return hi;
}

void PreProcess()
{
    bool convexA = isConvexPolygon(polygon1);
    bool convexB = isConvexPolygon(polygon2);
    if ((convexA && convexB) || (n1 <= 100 && n2 <= 100))
    {
        useNfpMode = true;
        g_nfpPoly = nfpCombined(polygon1, polygon2);
        if (!polygon2.vertices.empty())
            g_refB = polygon2.vertices[0];
        return;
    }

    useNfpMode = false;
    const int kDirectionCount = 8;
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> angleGen(0.0, 2.0 * PI);
    double baseTheta = angleGen(gen);
    g_escapeDirs = BuildUniformDirections(kDirectionCount, baseTheta);
    g_escapeDist.resize(kDirectionCount);

    for (int i = 0; i < kDirectionCount; ++i)
    {
        Vector2D unit = g_escapeDirs[i].Normalize();
        g_escapeDist[i] = GenSolution(unit);
        g_escapeDirs[i] = unit * g_escapeDist[i];
    }
}

// ------------------------
// 主函数
// ------------------------

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    std::string sok;

    if (!(std::cin >> n1 >> n2))
        return 1;
    if (n1 <= 0 || n2 <= 0)
    {
        cerr << "Input data error: the number of vertices of both polygons should be greater than 2" << endl;
        return 1;
    }

    polygon1.vertices.resize(n1);
    for (int i = 0; i < n1; ++i)
    {
        cin >> polygon1.vertices[i].x >> polygon1.vertices[i].y;
    }
    polygon2.vertices.resize(n2);
    for (int i = 0; i < n2; ++i)
    {
        cin >> polygon2.vertices[i].x >> polygon2.vertices[i].y;
    }

    cin >> sok;
    if (sok != "OK")
    {
        cerr << "Input data error: waiting for OK after obtaining polygons but I get " << sok << endl;
        return 0;
    }

    PreProcess();
    cout << "OK" << endl;
    cout.flush();

    cin >> m;
    testCases.resize(m);
    for (int i = 0; i < m; ++i)
    {
        cin >> testCases[i].x >> testCases[i].y;
    }

    cin >> sok;
    if (sok != "OK")
    {
        cerr << "Input data error: waiting for OK after that I have received all test points but I get " << sok << endl;
        return 0;
    }
    cout << m << '\n';

    if (useNfpMode)
    {
        vector<Vector2D> results;
        results.reserve(m);
        for (int i = 0; i < m; ++i)
        {
            Vector2D newRef = g_refB + testCases[i];
            DistanceResult distResult = pointToPolygon(newRef, g_nfpPoly);
            bool needsSeparation = distResult.inside && distResult.distance > EPS;
            if (needsSeparation)
                results.push_back(distResult.nearest - newRef);
            else
                results.emplace_back(0, 0);
        }
        for (int i = 0; i < m; ++i)
        {
            cout << fixed << setprecision(5) << results[i].x << " " << results[i].y << '\n';
            cout.flush();
        }
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            Vector2D res;
            double anslen = std::numeric_limits<double>::infinity();
            for (int j = 0; j < (int)g_escapeDirs.size(); ++j)
            {
                Vector2D tmp = g_escapeDirs[j] - testCases[i];
                double dist = tmp.Length();
                if (dist < anslen)
                {
                    anslen = dist;
                    res = tmp;
                }
            }
            cout << fixed << setprecision(5) << res.x << " " << res.y << '\n';
            cout.flush();
        }
    }

    cout << "OK" << '\n';
    cout.flush();
    return 0;
}
