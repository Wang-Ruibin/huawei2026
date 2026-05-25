#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <unistd.h>
#include <thread>
#include <chrono>

struct Vector2D;
struct Polygon;

using namespace std;
const double EPS = 1e-6;

struct Vector2D
{
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}
    Vector2D operator-(const Vector2D &other) const { return Vector2D(x - other.x, y - other.y); }
    Vector2D operator+(const Vector2D &other) const { return Vector2D(x + other.x, y + other.y); }
    Vector2D operator*(double scalar) const { return Vector2D(x * scalar, y * scalar); }
    double Dot(const Vector2D &other) const { return x * other.x + y * other.y; }
    double Length() const { return std::sqrt(x * x + y * y); }
    Vector2D Normalize() const
    {
        double len = Length();
        if (len == 0)
        {
            return *this;
        }
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

struct Polygon
{
    std::vector<Vector2D> vertices;

    Polygon() = default;
    Polygon(std::initializer_list<Vector2D> vts) : vertices(vts) {}
    Vector2D GetCenter() const
    {
        Vector2D center(0, 0);
        if (vertices.empty())
        {
            return center;
        }
        for (const auto &v : vertices)
        {
            center = center + v;
        }
        return center * (1.0 / vertices.size());
    }
    void MoveByVec(const Vector2D &vec)
    {
        for (auto &v : vertices)
        {
            v = v + vec;
        }
    }
};

double Cross(const Vector2D &a, const Vector2D &b)
{
    return a.x * b.y - a.y * b.x;
}

double Cross(const Vector2D &o, const Vector2D &a, const Vector2D &b)
{
    return Cross(a - o, b - o);
}

vector<Vector2D> ConvexHull(vector<Vector2D> pts)
{
    int n = pts.size();
    if (n <= 2)
        return pts;
    sort(pts.begin(), pts.end(), [](const Vector2D &a, const Vector2D &b)
         {
             if (fabs(a.x - b.x) > EPS) return a.x < b.x;
             return a.y < b.y; });

    vector<Vector2D> hull;

    for (int i = 0; i < n; ++i)
    {
        while (hull.size() >= 2 &&
               Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= 0)
        {
            hull.pop_back();
        }
        hull.push_back(pts[i]);
    }

    int t = hull.size();
    for (int i = n - 2; i >= 0; --i)
    {
        while ((int)hull.size() > t &&
               Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= 0)
        {
            hull.pop_back();
        }
        hull.push_back(pts[i]);
    }

    hull.pop_back();
    return hull;
}

int n1 = 0, n2 = 0, m = 0;
bool f1, f2;
Polygon polygon1;
Polygon polygon2;
Polygon polygon3, polygon4;
vector<Vector2D> testCases;
vector<Vector2D> axes;

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
        {
            minProj = proj;
        }
        if (proj > maxProj)
        {
            maxProj = proj;
        }
    }
    return {minProj, maxProj};
}

double calOverlapWithDir(const Projection &projA, const Projection &projB, const Vector2D &axis, Vector2D &outAxis)
{
    // 分离
    if (projA.max <= projB.min || projB.max <= projA.min)
        return 0.0;

    double d1 = projA.max - projB.min;
    double d2 = projB.max - projA.min;

    if (d1 < d2)
    {
        outAxis = axis; // 保持方向
        return d1;
    }
    else
    {
        outAxis = axis * -1.0; // 反方向
        return d2;
    }
}

Vector2D GenSolution(const Vector2D &vec)
{
    Polygon polyB = polygon2;
    polyB.MoveByVec(vec);

    double minOverlap = std::numeric_limits<double>::infinity();
    Vector2D smallestAxis, outAxis;

    for (const auto &axis : axes)
    {
        Projection projA, projB;
        projA.min = polygon1.vertices[0].Dot(axis);
        projA.max = polygon1.vertices[1].Dot(axis);
        if (projA.min > projA.max)
            std::swap(projA.min, projA.max);
        projB.min = polyB.vertices[0].Dot(axis);
        projB.max = polyB.vertices[1].Dot(axis);
        if (projB.min > projB.max)
            std::swap(projB.min, projB.max);
        double overlap = calOverlapWithDir(projA, projB, axis, outAxis);

        for (size_t i = 2; i < std::max(polygon1.vertices.size(), polyB.vertices.size()); ++i)
        {
            if (overlap >= minOverlap)
                break;
            if (i < polygon1.vertices.size())
            {
                double proj = polygon1.vertices[i].Dot(axis);
                if (proj < projA.min)
                    projA.min = proj;
                if (proj > projA.max)
                    projA.max = proj;
            }
            if (i < polyB.vertices.size())
            {
                double proj = polyB.vertices[i].Dot(axis);
                if (proj < projB.min)
                    projB.min = proj;
                if (proj > projB.max)
                    projB.max = proj;
            }
            overlap = calOverlapWithDir(projA, projB, axis, outAxis);
        }

        if (projA.max <= projB.min || projB.max <= projA.min)
        {
            return {0.0, 0.0};
        }

        if (overlap < minOverlap)
        {
            minOverlap = overlap;
            smallestAxis = outAxis;
        }
    }

    return {smallestAxis * minOverlap};
}

// 选手在规定的时间内进行预处理，完成后返回OK
void PreProcess()
{
    polygon3.vertices = ConvexHull(polygon1.vertices);
    polygon4.vertices = ConvexHull(polygon2.vertices);

    auto addAxes = [&](const Polygon &poly)
    {
        int n = poly.vertices.size();
        for (int i = 0; i < n; ++i)
        {
            Vector2D p1 = poly.vertices[i];
            Vector2D p2 = poly.vertices[(i + 1) % n];
            Vector2D edge = p2 - p1;

            Vector2D axis = edge.Perp().Normalize();

            if (axis.x < 0 || (fabs(axis.x) < EPS && axis.y < 0))
                axis = axis * -1.0;

            axes.push_back(axis);
        }
    };

    addAxes(polygon1);
    addAxes(polygon2);

    if (polygon3.vertices.size() != polygon1.vertices.size())
        addAxes(polygon3), f1 = true;
    if (polygon4.vertices.size() != polygon2.vertices.size())
        addAxes(polygon4), f2 = true;

    sort(axes.begin(), axes.end(), [](const Vector2D &a, const Vector2D &b)
         { return atan2(a.y, a.x) < atan2(b.y, b.x); });

    vector<Vector2D> uniqueAxes;
    for (const auto &axis : axes)
    {
        if (uniqueAxes.empty())
        {
            uniqueAxes.push_back(axis);
        }
        else
        {
            Vector2D last = uniqueAxes.back();
            if (fabs(last.x * axis.y - last.y * axis.x) > EPS)
            {
                uniqueAxes.push_back(axis);
            }
        }
    }
    axes = uniqueAxes;
}

int main()
{
    // freopen("practice_3.in", "r", stdin);
    // freopen("practice_3.out", "w", stdout);
    // =============== 1. read polygons ===================
    cin >> n1 >> n2;

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

    // wait for OK to ensure all polygon data is received
    string okResp;
    cin >> okResp;
    if (okResp != "OK")
    {
        cerr << "Input data error: waiting for OK after obtaining polygons but I get " << okResp << endl;
        return 0;
    }

    // ============== 2. preprocess ===================
    PreProcess();
    // cout << "f1 " << f1 << " f2 " << f2 << endl;
    //  send OK after finishing preprocessing
    cout << "OK" << endl;
    cout.flush();

    // ============== 3. read test data ===================
    cin >> m;
    testCases.resize(m);

    for (int i = 0; i < m; ++i)
    {
        cin >> testCases[i].x >> testCases[i].y;
    }

    // wait for OK to ensure all test cases are received
    cin >> okResp;
    if (okResp != "OK")
    {
        cerr << "Input data error: waiting for OK after that I have received all test points but I get " << okResp
             << endl;
        return 0;
    }

    // std::this_thread::sleep_for(std::chrono::seconds(15));
    //  ================ 4. solve and output results ===================

    cout << m << "\n";
    for (int i = 0; i < m; ++i)
    {
        const Vector2D &res = GenSolution(testCases[i]);
        cout << fixed << std::setprecision(5) << res.x << " " << res.y << "\n";
        cout.flush();
    }

    // send OK after outputting all answer
    cout << "OK" << endl;
    cout.flush();

    return 0;
}

// int main()
// {
//     freopen("practice_1.in", "r", stdin);
//     freopen("practice_1.out", "w", stdout);
//     // =============== 1. read polygons ===================
//     cin >> n1 >> n2;

//     if (n1 <= 0 || n2 <= 0)
//     {
//         cerr << "Input data error: the number of vertices of both polygons should be greater than 2" << endl;
//         return 1;
//     }

//     polygon1.vertices.resize(n1);
//     for (int i = 0; i < n1; ++i)
//     {
//         cin >> polygon1.vertices[i].x >> polygon1.vertices[i].y;
//     }

//     polygon2.vertices.resize(n2);
//     for (int i = 0; i < n2; ++i)
//     {
//         cin >> polygon2.vertices[i].x >> polygon2.vertices[i].y;
//     }

//     // wait for OK to ensure all polygon data is received
//     // string okResp;
//     // cin >> okResp;
//     // if (okResp != "OK")
//     // {
//     //     cerr << "Input data error: waiting for OK after obtaining polygons but I get " << okResp << endl;
//     //     return 0;
//     // }

//     // ============== 2. preprocess ===================
//     PreProcess();
//     // send OK after finishing preprocessing
//     // cout << "OK" << endl;
//     // cout.flush();

//     // ============== 3. read test data ===================
//     cin >> m;
//     testCases.resize(m);

//     for (int i = 0; i < m; ++i)
//     {
//         cin >> testCases[i].x >> testCases[i].y;
//     }

//     // wait for OK to ensure all test cases are received
//     // cin >> okResp;
//     // if (okResp != "OK")
//     // {
//     //     cerr << "Input data error: waiting for OK after that I have received all test points but I get " << okResp
//     //          << endl;
//     //     return 0;
//     // }

//     // ================ 4. solve and output results ===================
//     for (int i = 0; i < m; ++i)
//     {
//         const Vector2D &res = GenSolution(testCases[i]);
//         cout << fixed << std::setprecision(5) << res.x << " " << res.y << endl;
//         cout.flush();
//     }

//     // send OK after outputting all answer
//     // cout << "OK" << endl;
//     // cout.flush();

//     return 0;
// }