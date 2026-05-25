#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace std;

const double EPS = 1e-9;

struct Vector2D
{
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}

    Vector2D operator-(const Vector2D &o) const { return {x - o.x, y - o.y}; }
    Vector2D operator+(const Vector2D &o) const { return {x + o.x, y + o.y}; }
    Vector2D operator*(double s) const { return {x * s, y * s}; }

    double Dot(const Vector2D &o) const { return x * o.x + y * o.y; }
    double Length() const { return sqrt(x * x + y * y); }

    Vector2D Normalize() const
    {
        double len = Length();
        if (len < EPS)
            return {0, 0};
        return {x / len, y / len};
    }

    Vector2D Perp() const { return {-y, x}; }
};

struct Polygon
{
    vector<Vector2D> vertices;

    Vector2D GetCenter() const
    {
        Vector2D c(0, 0);
        for (auto &v : vertices)
            c = c + v;
        return c * (1.0 / vertices.size());
    }

    void MoveByVec(const Vector2D &v)
    {
        for (auto &p : vertices)
            p = p + v;
    }
};

struct Projection
{
    double min, max;
};

struct Piece
{
    Polygon poly;
    vector<Vector2D> axes;
    Vector2D center;
};

Polygon polygon1, polygon2;
vector<Piece> parts1, parts2;

double Cross(const Vector2D &a, const Vector2D &b)
{
    return a.x * b.y - a.y * b.x;
}

double Cross(const Vector2D &o, const Vector2D &a, const Vector2D &b)
{
    return Cross(a - o, b - o);
}

double SignedArea(const vector<Vector2D> &p)
{
    double a = 0;
    for (int i = 0; i < p.size(); i++)
        a += Cross(p[i], p[(i + 1) % p.size()]);
    return a * 0.5;
}

void EnsureCCW(vector<Vector2D> &p)
{
    if (SignedArea(p) < 0)
        reverse(p.begin(), p.end());
}

bool PointInTriangle(const Vector2D &p, const Vector2D &a, const Vector2D &b, const Vector2D &c)
{
    double c1 = Cross(a, b, p);
    double c2 = Cross(b, c, p);
    double c3 = Cross(c, a, p);
    return c1 >= -EPS && c2 >= -EPS && c3 >= -EPS;
}

bool IsEar(int i, const vector<Vector2D> &poly)
{
    int n = poly.size();
    int pre = (i - 1 + n) % n;
    int nxt = (i + 1) % n;

    if (Cross(poly[pre], poly[i], poly[nxt]) <= EPS)
        return false;

    for (int j = 0; j < n; j++)
    {
        if (j == pre || j == i || j == nxt)
            continue;
        if (PointInTriangle(poly[j], poly[pre], poly[i], poly[nxt]))
            return false;
    }

    return true;
}

vector<Polygon> EarClip(const Polygon &poly)
{
    vector<Polygon> res;
    vector<Vector2D> pts = poly.vertices;

    EnsureCCW(pts);

    while (pts.size() > 3)
    {
        bool found = false;
        for (int i = 0; i < pts.size(); i++)
        {
            if (IsEar(i, pts))
            {
                int pre = (i - 1 + pts.size()) % pts.size();
                int nxt = (i + 1) % pts.size();

                res.push_back({{pts[pre], pts[i], pts[nxt]}});
                pts.erase(pts.begin() + i);
                found = true;
                break;
            }
        }
        if (!found)
            break;
    }

    if (pts.size() == 3)
        res.push_back({{pts[0], pts[1], pts[2]}});

    return res;
}

vector<Vector2D> CollectAxes(const Polygon &poly)
{
    vector<Vector2D> axes;

    for (int i = 0; i < poly.vertices.size(); i++)
    {
        auto e = poly.vertices[(i + 1) % poly.vertices.size()] - poly.vertices[i];
        auto axis = e.Perp().Normalize();

        if (axis.x < 0 || (fabs(axis.x) < EPS && axis.y < 0))
            axis = axis * -1;

        axes.push_back(axis);
    }

    return axes;
}

Projection Project(const Polygon &poly, const Vector2D &axis)
{
    double mn = poly.vertices[0].Dot(axis);
    double mx = mn;

    for (int i = 1; i < poly.vertices.size(); i++)
    {
        double v = poly.vertices[i].Dot(axis);
        mn = min(mn, v);
        mx = max(mx, v);
    }

    return {mn, mx};
}

// ⭐⭐⭐ 核心：最小分离距离（支持包含）
double PenetrationDepth(const Projection &A, const Projection &B)
{
    if (A.max <= B.min + EPS || B.max <= A.min + EPS)
        return 0.0;

    double pushLeft = A.max - B.min;
    double pushRight = B.max - A.min;

    return min(pushLeft, pushRight);
}

bool SAT(const Polygon &A, const Polygon &B, double &depth, Vector2D &axis)
{
    vector<Vector2D> axes = CollectAxes(A);
    auto bAxes = CollectAxes(B);
    axes.insert(axes.end(), bAxes.begin(), bAxes.end());

    depth = numeric_limits<double>::infinity();

    for (auto &ax : axes)
    {
        auto pa = Project(A, ax);
        auto pb = Project(B, ax);

        if (pa.max <= pb.min + EPS || pb.max <= pa.min + EPS)
            return false;

        double d = PenetrationDepth(pa, pb);

        if (d < depth)
        {
            depth = d;
            axis = ax;
        }
    }

    return true;
}

Vector2D GenSolution(const Vector2D &vec)
{
    Polygon B = polygon2;
    B.MoveByVec(vec);

    auto movedParts2 = parts2;
    for (auto &p : movedParts2)
    {
        p.poly.MoveByVec(vec);
        p.center = p.poly.GetCenter();
    }

    double best = numeric_limits<double>::infinity();
    Vector2D bestAxis;

    bool hit = false;

    for (auto &A : parts1)
    {
        for (auto &B : movedParts2)
        {
            double depth;
            Vector2D axis;

            if (!SAT(A.poly, B.poly, depth, axis))
                continue;

            hit = true;

            Vector2D dir = B.center - A.center;
            if (axis.Dot(dir) < 0)
                axis = axis * -1;

            if (depth < best)
            {
                best = depth;
                bestAxis = axis;
            }
        }
    }

    if (!hit)
        return {0, 0};

    return bestAxis * best;
}

void PreProcess()
{
    auto tris1 = EarClip(polygon1);
    auto tris2 = EarClip(polygon2);

    for (auto &t : tris1)
        parts1.push_back({t, CollectAxes(t), t.GetCenter()});

    for (auto &t : tris2)
        parts2.push_back({t, CollectAxes(t), t.GetCenter()});
}

int main()
{
    int n1, n2, m;
    cin >> n1 >> n2;

    polygon1.vertices.resize(n1);
    polygon2.vertices.resize(n2);

    for (auto &v : polygon1.vertices)
        cin >> v.x >> v.y;
    for (auto &v : polygon2.vertices)
        cin >> v.x >> v.y;

    string ok;
    cin >> ok;

    PreProcess();

    cout << "OK\n";

    cin >> m;
    vector<Vector2D> tests(m);
    for (auto &t : tests)
        cin >> t.x >> t.y;

    cin >> ok;

    for (auto &t : tests)
    {
        auto r = GenSolution(t);
        cout << fixed << setprecision(5) << r.x << " " << r.y << "\n";
    }

    cout << "OK\n";
}