#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <string>
#include <limits>

using namespace std;

constexpr double EPS = 1e-9;
constexpr double EPS2 = EPS * EPS;
constexpr double INF = 1e20;

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

inline double pointSegDist2(const Pt &p, const Pt &a, const Pt &b, Pt &q)
{
    Vec ab = b - a;
    Vec ap = p - a;
    double tDen = dot(ab, ab);

    if (tDen < EPS2)
    {
        q = a;
        return dot(ap, ap);
    }

    double t = dot(ap, ab) / tDen;
    if (t <= 0.0)
    {
        q = a;
        return dot(ap, ap);
    }
    if (t >= 1.0)
    {
        q = b;
        Vec bp = p - b;
        return dot(bp, bp);
    }

    q = a + ab * t;
    Vec pq = p - q;
    return dot(pq, pq);
}

inline double pseudoAngle(double x, double y)
{
    const double s = 1e-12;
    if (absd(x) < s && absd(y) < s)
        return 0.0;
    if (x >= 0 && y >= 0)
        return (absd(x) < s) ? 1.0 : y / (x + y + s);
    if (x < 0 && y >= 0)
        return 1.0 + (absd(y) < s ? 1.0 : -x / (-x + y + s));
    if (x < 0 && y < 0)
        return 2.0 + (absd(x) < s ? 1.0 : -y / (-x - y + s));
    return 3.0 + (absd(y) < s ? 1.0 : x / (x - y + s));
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

struct Closest
{
    Pt q;
    double d;
    bool inside;
};

// ===================== BVH for closest-edge query + point-in-poly =====================

struct BBox
{
    double minX, maxX, minY, maxY;
};

struct BVHEdge
{
    Pt a, b;
    BBox box;
};

struct BVHNode
{
    BBox box;
    int left = -1, right = -1;
    int l = 0, r = 0;
    bool leaf = false;
};

static vector<BVHEdge> g_edges;
static vector<int> g_order;
static vector<BVHNode> g_tree;
static const Polygon *g_cachedPoly = nullptr;
static int g_root = -1;

inline BBox mergeBox(const BBox &A, const BBox &B)
{
    return {
        min(A.minX, B.minX),
        max(A.maxX, B.maxX),
        min(A.minY, B.minY),
        max(A.maxY, B.maxY)};
}

inline BBox boxOfSeg(const Pt &a, const Pt &b)
{
    return {
        min(a.x, b.x),
        max(a.x, b.x),
        min(a.y, b.y),
        max(a.y, b.y)};
}

inline double boxDist2(const Pt &p, const BBox &b)
{
    double dx = 0.0, dy = 0.0;

    if (p.x < b.minX)
        dx = b.minX - p.x;
    else if (p.x > b.maxX)
        dx = p.x - b.maxX;

    if (p.y < b.minY)
        dy = b.minY - p.y;
    else if (p.y > b.maxY)
        dy = p.y - b.maxY;

    return dx * dx + dy * dy;
}

inline bool isConvexPoly(const Polygon &g)
{
    int n = g.n;
    if (n < 3)
        return false;

    int sgn = 0;
    for (int i = 0; i < n; ++i)
    {
        const Pt &a = g.v[i];
        const Pt &b = g.v[(i + 1) % n];
        const Pt &c = g.v[(i + 2) % n];
        double z = cross(b - a, c - b);
        if (absd(z) < EPS)
            continue;
        int cur = (z > 0 ? 1 : -1);
        if (sgn == 0)
            sgn = cur;
        else if (sgn != cur)
            return false;
    }
    return true;
}

inline bool pointInConvexPoly(const vector<Pt> &p, const Pt &q)
{
    int n = (int)p.size();
    if (n < 3)
        return false;

    if (cross(p[1] - p[0], q - p[0]) < -EPS)
        return false;
    if (cross(p[n - 1] - p[0], q - p[0]) > EPS)
        return false;

    int l = 1, r = n - 1;
    while (r - l > 1)
    {
        int mid = (l + r) >> 1;
        if (cross(p[mid] - p[0], q - p[0]) >= 0)
            l = mid;
        else
            r = mid;
    }

    return cross(p[(l + 1) % n] - p[l], q - p[l]) >= -EPS;
}

Closest closestOnConvexPolygon(const Pt &p, const Polygon &g)
{
    Closest ans;
    ans.inside = pointInConvexPoly(g.v, p);
    ans.q = p;
    ans.d = INF;

    if (!ans.inside || g.n < 3)
        return ans;

    double best = INF;
    Pt bestQ = p;
    double bestY = -INF, bestX = -INF;

    for (int i = 0; i < g.n; ++i)
    {
        Pt q;
        double d2 = pointSegDist2(p, g.v[i], g.v[(i + 1) % g.n], q);

        Vec u = q - p;
        double uy = u.y, ux = u.x;

        if (d2 < best - EPS2)
        {
            best = d2;
            bestQ = q;
            bestY = uy;
            bestX = ux;
        }
        else if (absd(d2 - best) < EPS2)
        {
            if (uy > bestY + EPS || (absd(uy - bestY) < EPS && ux > bestX + EPS))
            {
                bestQ = q;
                bestY = uy;
                bestX = ux;
            }
        }
    }

    ans.q = bestQ;
    ans.d = sqrt(best);
    return ans;
}

int buildBVH(int l, int r)
{
    int id = (int)g_tree.size();
    g_tree.push_back(BVHNode());

    BBox box = g_edges[g_order[l]].box;
    for (int i = l + 1; i < r; ++i)
        box = mergeBox(box, g_edges[g_order[i]].box);

    g_tree[id].box = box;
    g_tree[id].l = l;
    g_tree[id].r = r;

    if (r - l <= 8)
    {
        g_tree[id].leaf = true;
        return id;
    }

    double spanX = box.maxX - box.minX;
    double spanY = box.maxY - box.minY;
    int axis = (spanX >= spanY ? 0 : 1);
    int mid = (l + r) >> 1;

    nth_element(g_order.begin() + l, g_order.begin() + mid, g_order.begin() + r,
                [&](int i, int j)
                {
                    const BVHEdge &ei = g_edges[i];
                    const BVHEdge &ej = g_edges[j];
                    double ci = axis == 0 ? (ei.a.x + ei.b.x) * 0.5 : (ei.a.y + ei.b.y) * 0.5;
                    double cj = axis == 0 ? (ej.a.x + ej.b.x) * 0.5 : (ej.a.y + ej.b.y) * 0.5;
                    return ci < cj;
                });

    g_tree[id].left = buildBVH(l, mid);
    g_tree[id].right = buildBVH(mid, r);
    return id;
}

void BuildClosestBVH(const Polygon &poly)
{
    if (g_cachedPoly == &poly)
        return;

    g_cachedPoly = &poly;
    g_edges.clear();
    g_order.clear();
    g_tree.clear();
    g_root = -1;

    if (poly.n < 3)
        return;

    g_edges.reserve(poly.n);
    g_order.reserve(poly.n);

    for (int i = 0; i < poly.n; ++i)
    {
        Pt a = poly.v[i];
        Pt b = poly.v[(i + 1) % poly.n];
        g_edges.push_back({a, b, boxOfSeg(a, b)});
        g_order.push_back(i);
    }

    g_tree.reserve(poly.n * 2);
    g_root = buildBVH(0, (int)g_order.size());
}

inline void relaxClosest(const Pt &p, const Pt &q, double d2, Closest &ans,
                         double &best, double &bestY, double &bestX)
{
    Vec u = q - p;
    double uy = u.y;
    double ux = u.x;

    if (d2 < best - EPS2)
    {
        best = d2;
        ans.q = q;
        bestY = uy;
        bestX = ux;
    }
    else if (absd(d2 - best) < EPS2)
    {
        if (uy > bestY + EPS || (absd(uy - bestY) < EPS && ux > bestX + EPS))
        {
            ans.q = q;
            bestY = uy;
            bestX = ux;
        }
    }
}

void queryBVH(int id, const Pt &p, Closest &ans, double &best, double &bestY, double &bestX)
{
    const BVHNode &node = g_tree[id];

    if (boxDist2(p, node.box) > best + EPS2)
        return;

    if (node.leaf)
    {
        for (int i = node.l; i < node.r; ++i)
        {
            const BVHEdge &e = g_edges[g_order[i]];
            Pt q;
            double d2 = pointSegDist2(p, e.a, e.b, q);
            relaxClosest(p, q, d2, ans, best, bestY, bestX);
        }
        return;
    }

    int L = node.left, R = node.right;
    double dL = boxDist2(p, g_tree[L].box);
    double dR = boxDist2(p, g_tree[R].box);

    if (dL < dR)
    {
        if (dL <= best + EPS2)
            queryBVH(L, p, ans, best, bestY, bestX);
        if (dR <= best + EPS2)
            queryBVH(R, p, ans, best, bestY, bestX);
    }
    else
    {
        if (dR <= best + EPS2)
            queryBVH(R, p, ans, best, bestY, bestX);
        if (dL <= best + EPS2)
            queryBVH(L, p, ans, best, bestY, bestX);
    }
}

bool pointInPolyBVHRec(int id, const Pt &p, bool &onEdge)
{
    const BVHNode &node = g_tree[id];

    if (p.y < node.box.minY - EPS || p.y > node.box.maxY + EPS || p.x > node.box.maxX + EPS)
        return false;

    if (node.leaf)
    {
        bool inside = false;
        for (int i = node.l; i < node.r; ++i)
        {
            const BVHEdge &e = g_edges[g_order[i]];
            const Pt &a = e.a;
            const Pt &b = e.b;

            if (pointOnSeg(p, a, b))
            {
                onEdge = true;
                return true;
            }

            if ((a.y > p.y) != (b.y > p.y))
            {
                double x = a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y);
                if (x > p.x + EPS)
                    inside = !inside;
            }
        }
        return inside;
    }

    bool left = pointInPolyBVHRec(node.left, p, onEdge);
    if (onEdge)
        return true;
    bool right = pointInPolyBVHRec(node.right, p, onEdge);
    return left ^ right;
}

bool pointInPolyBVH(const Pt &p, const Polygon &g)
{
    BuildClosestBVH(g);
    if (g_root == -1)
        return false;

    const BBox &box = g_tree[g_root].box;
    if (p.x < box.minX - EPS || p.x > box.maxX + EPS || p.y < box.minY - EPS || p.y > box.maxY + EPS)
        return false;

    bool onEdge = false;
    bool inside = pointInPolyBVHRec(g_root, p, onEdge);
    return inside || onEdge;
}

Closest closestOnPolygon(const Pt &p, const Polygon &g)
{
    BuildClosestBVH(g);

    Closest ans;
    ans.inside = pointInPolyBVH(p, g);
    ans.d = INF;
    ans.q = p;

    if (g_root == -1 || !ans.inside)
        return ans;

    double best = INF;
    double bestY = -INF;
    double bestX = -INF;

    queryBVH(g_root, p, ans, best, bestY, bestX);

    if (best < INF)
    {
        ans.d = sqrt(best);
        if (best < EPS2)
            ans.q = p;
    }

    return ans;
}

// ===================== NFP =====================

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
            ea.push_back({e, pseudoAngle(e.x, e.y)});
    }

    for (int i = 0; i < b.n; ++i)
    {
        Vec e = b.v[(sb + i + 1) % b.n] - b.v[(sb + i) % b.n];
        e = -e;
        if (len2(e) > EPS2)
            eb.push_back({e, pseudoAngle(e.x, e.y)});
    }

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
    buffer.reserve(m * 32);

    bool nfpConv = isConvexPoly(nfp);

    Pt ref = polygon2.v.empty() ? Pt(0, 0) : polygon2.v[0];
    for (int i = 0; i < m; ++i)
    {
        const Pt &t = testCases[i];
        Pt cur(ref.x + t.x, ref.y + t.y);

        Closest ans;
        if (nfpConv)
            ans = closestOnConvexPolygon(cur, nfp);
        else
            ans = closestOnPolygon(cur, nfp); // 你原来的 BVH 版本

        if (ans.inside && ans.d > EPS)
        {
            Pt out = ans.q - cur;
            char line[64];
            snprintf(line, sizeof(line), "%.5f %.5f\n", out.x, out.y);
            buffer += line;
        }
        else
        {
            buffer += "0.00000 0.00000\n";
        }
    }
    cout << buffer;
    cout << "OK\n";
    cout.flush();
    return 0;
}