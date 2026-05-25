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
#include <random>

static std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

struct Vector2D;
struct Polygon;

using namespace std;
const double EPS = 1e-8; // 些许提升，可以再看看
const double eps = 1e-9;

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
vector<Vector2D> axes, allaxes;

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

double mind1, mind2; // d1为正方向，d2要转向

double calOverlapWithDir(const Projection &projA, const Projection &projB)
{
    // 分离
    if (projA.max <= projB.min || projB.max <= projA.min)
        return 0.0;

    double d1 = projA.max - projB.min;
    double d2 = projB.max - projA.min;

    if (d1 < d2)
    {
        return d1;
    }
    else
    {
        return d2;
    }
}

bool OnSegment(const Vector2D &A, const Vector2D &B, const Vector2D &P)
{
    // 1. 检查共线性：叉积必须为0（允许浮点误差）
    if (std::fabs(Cross(A, B, P)) > eps)
        return false;

    // 2. 检查投影范围：点P在线段AB的包围盒内
    //    使用坐标范围判断，避免点积可能出现的精度问题
    double minX = std::min(A.x, B.x) - eps;
    double maxX = std::max(A.x, B.x) + eps;
    double minY = std::min(A.y, B.y) - eps;
    double maxY = std::max(A.y, B.y) + eps;

    return (P.x >= minX && P.x <= maxX && P.y >= minY && P.y <= maxY);
}

bool PointInPolygonStrict(const Vector2D &p, const Polygon &poly)
{
    const auto &v = poly.vertices;
    int n = v.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        const Vector2D &a = v[i];
        const Vector2D &b = v[j];
        // 点在边上 → 不在内部
        if (OnSegment(a, b, p))
            return false;
        // 射线法：从p向右发射水平射线
        if ((a.y > p.y) != (b.y > p.y))
        {
            double x = a.x + (b.x - a.x) * (p.y - a.y) / (b.y - a.y);
            if (x > p.x)
                inside = !inside;
        }
    }
    return inside;
}

int SegmentsStrictIntersect(const Vector2D &a, const Vector2D &b,
                            const Vector2D &c, const Vector2D &d)
{
    double d1 = Cross(a, b, c);
    double d2 = Cross(a, b, d);
    double d3 = Cross(c, d, a);
    double d4 = Cross(c, d, b);

    // 规范相交（交点在线段内部）
    if (d1 * d2 < -EPS && d3 * d4 < -EPS)
        return 5;

    // 端点相交
    if (OnSegment(c, d, a))
        return 1;
    if (OnSegment(c, d, b))
        return 2;
    if (OnSegment(a, b, c))
        return 3;
    if (OnSegment(a, b, d))
        return 4;

    return 0;
}

struct edgeacross
{
    Vector2D a, b, c, d;
};

edgeacross edge;
vector<Vector2D> point;

// 获取多边形的包围盒（最小/最大坐标）
struct BBox
{
    double minX, maxX, minY, maxY;
};

BBox GetBBox(const Polygon &poly)
{
    if (poly.vertices.empty())
    {
        return {0, 0, 0, 0};
    }
    double minX = poly.vertices[0].x, maxX = poly.vertices[0].x;
    double minY = poly.vertices[0].y, maxY = poly.vertices[0].y;
    for (const auto &v : poly.vertices)
    {
        if (v.x < minX)
            minX = v.x;
        if (v.x > maxX)
            maxX = v.x;
        if (v.y < minY)
            minY = v.y;
        if (v.y > maxY)
            maxY = v.y;
    }
    return {minX, maxX, minY, maxY};
}

bool bboxOverlap(const Polygon &P, const Polygon &Q)
{
    BBox bboxP = GetBBox(P);
    BBox bboxQ = GetBBox(Q);
    // 使用 eps 容差，允许边界恰好接触的情况
    return !(bboxP.maxX + eps < bboxQ.minX || bboxQ.maxX + eps < bboxP.minX ||
             bboxP.maxY + eps < bboxQ.minY || bboxQ.maxY + eps < bboxP.minY);
}

Vector2D movepoint(Vector2D node, Vector2D vec)
{
    Vector2D tmp = node + vec;
    return tmp;
}

bool check(const Polygon &P, const Polygon &Q, const Vector2D &vec)
{
    // 1. 包围盒快速排除
    if (!bboxOverlap(P, Q))
        return true; // 需实现bboxOverlap

    // 2. 检查边严格相交
    // for (auto e : edge)
    // if (SegmentsStrictIntersect(edge.a, edge.b, edge.c + vec, edge.d + vec))
    //     return false;
    // edge.clear();

    const auto &vP = P.vertices;
    const auto &vQ = Q.vertices;
    int nP = vP.size(), nQ = vQ.size();
    for (int i = 0; i < nP; ++i)
    {
        const Vector2D &a = vP[i];
        const Vector2D &b = vP[(i + 1) % nP];
        for (int j = 0; j < nQ; ++j)
        {
            const Vector2D &c = vQ[j];
            const Vector2D &d = vQ[(j + 1) % nQ];
            Vector2D vec, point;
            int op = SegmentsStrictIntersect(a, b, c, d);
            if (op)
                return false;
            // if (op == 5)
            // {
            //     //edge = (edgeacross){a, b, c - vec, d - vec};
            //     return false;
            //     //  存在严格内部交点 → 重叠
            // }
            // if (op == 1)
            // {
            //     vec = b - a;
            //     // vec = vec.Normalize();
            //     vec = vec * EPS;
            //     point = movepoint(a, vec);
            //     if (PointInPolygonStrict(point, Q))
            //         return false;
            // }
            // if (op == 2)
            // {
            //     vec = a - b;
            //     // vec = vec.Normalize();
            //     vec = vec * EPS;
            //     point = movepoint(b, vec);
            //     if (PointInPolygonStrict(point, Q))
            //         return false;
            // }
            // if (op == 3)
            // {
            //     vec = d - c;
            //     // vec = vec.Normalize();
            //     vec = vec * EPS;
            //     point = movepoint(c, vec);
            //     if (PointInPolygonStrict(point, P))
            //         return false;
            // }
            // if (op == 4)
            // {
            //     vec = c - d;
            //     // vec = vec.Normalize();
            //     vec = vec * EPS;
            //     point = movepoint(d, vec);
            //     if (PointInPolygonStrict(point, P))
            //         return false;
            // }
        }
    }
    // if (!edge.empty())
    //     return false;
    // 3. 检查包含关系（一个多边形完全包含另一个，且无边接触）
    // 取P的任一顶点，若在Q内部则P包含于Q或相交，由于无边相交，只能是包含
    // for (auto node : point)
    //     if (PointInPolygonStrict(node + vec, P))
    //         return false;
    // point.clear();
    for (const auto &q : vQ)
    {
        if (PointInPolygonStrict(q, P))
            // point.push_back(q - vec);
            return false;
    }
    for (const auto &p : vP)
    {
        if (PointInPolygonStrict(p, Q))
            // point.push_back(q - vec);
            return false;
    }
    // if (!point.empty())
    //     return false;
    // 随机选点（顶点）进行包含检测
    // if (!vP.empty())
    // {
    //     int idx = std::uniform_int_distribution<int>(0, vP.size() - 1)(rng);
    //     if (PointInPolygonStrict(vP[idx], Q))
    //         return false;
    // }
    // if (!vQ.empty())
    // {
    //     int idx = std::uniform_int_distribution<int>(0, vQ.size() - 1)(rng);
    //     if (PointInPolygonStrict(vQ[idx], P))
    //         return false;
    // }

    // 4. 无重叠（分离或仅边界接触）
    return true;
}

Vector2D GenSolution(const Vector2D &Vec)
{
    Polygon polyB = polygon2, polyBT = polygon4;
    polyB.MoveByVec(Vec);
    polyBT.MoveByVec(Vec);

    if(check(polygon1, polyB, Vec))
        return {0.0, 0.0};

    double minOverlap = std::numeric_limits<double>::infinity();
    Vector2D smallestAxis;

    // for (auto &axis : axes)
    // {
    //     Projection projA, projB;
    //     projA.min = polygon3.vertices[0].Dot(axis);
    //     projA.max = polygon3.vertices[1].Dot(axis);
    //     if (projA.min > projA.max)
    //         std::swap(projA.min, projA.max);
    //     projB.min = polyBT.vertices[0].Dot(axis);
    //     projB.max = polyBT.vertices[1].Dot(axis);
    //     if (projB.min > projB.max)
    //         std::swap(projB.min, projB.max);
    //     double Overlap = calOverlapWithDir(projA, projB);

    //     for (size_t i = 2; i < std::max(polygon3.vertices.size(), polyBT.vertices.size()); ++i)
    //     {
    //         if (Overlap >= minOverlap)
    //             break;
    //         if (i < polygon3.vertices.size())
    //         {
    //             double proj = polygon3.vertices[i].Dot(axis);
    //             if (proj < projA.min)
    //                 projA.min = proj;
    //             if (proj > projA.max)
    //                 projA.max = proj;
    //         }
    //         if (i < polyBT.vertices.size())
    //         {
    //             double proj = polyBT.vertices[i].Dot(axis);
    //             if (proj < projB.min)
    //                 projB.min = proj;
    //             if (proj > projB.max)
    //                 projB.max = proj;
    //         }
    //         Overlap = calOverlapWithDir(projA, projB);
    //     }
    //     if (projA.max <= projB.min || projB.max <= projA.min)
    //         return {0.0, 0.0};
    //     if (Overlap < minOverlap)
    //     {
    //         minOverlap = Overlap;
    //         if (projA.max - projB.min <= projB.max - projA.min)
    //             smallestAxis = axis;
    //         else
    //             smallestAxis = axis * -1.0;
    //     }
    // }

    // if (f1 == false && f2 == false)
    //     return smallestAxis * minOverlap;

    for (const auto &axis : axes)
    {
        Polygon polyans1 = polyB, polyans2 = polyB;
        Vector2D vec = axis * minOverlap;
        polyans1.MoveByVec(vec);
        if (check(polygon1, polyans1, vec))
        {
            double maxn = minOverlap, minn = 0;
            double ans1 = maxn;
            while (maxn - minn > EPS)
            {
                double mid = (maxn + minn) / 2.0;
                Vector2D vect = axis * mid;
                Polygon polyans = polyB;
                polyans.MoveByVec(vect);
                if (check(polygon1, polyans, vect))
                    maxn = mid, ans1 = mid;
                else
                    minn = mid;
            }
            if (ans1 < minOverlap)
            {
                minOverlap = ans1;
                smallestAxis = axis;
            }
        }
        vec = axis * minOverlap;
        vec = vec * -1.0;
        polyans2.MoveByVec(vec);
        if (check(polygon1, polyans2, vec))
        {
            double maxn = minOverlap, minn = 0;
            double ans2 = maxn;
            while (maxn - minn > EPS)
            {
                double mid = (maxn + minn) / 2.0;
                Vector2D vect = axis * mid;
                vect = vect * -1.0;
                Polygon polyans = polyB;
                polyans.MoveByVec(vect);
                if (check(polygon1, polyans, vect))
                    maxn = mid, ans2 = mid;
                else
                    minn = mid;
            }
            if (ans2 < minOverlap)
            {
                minOverlap = ans2;
                smallestAxis = axis * -1.0;
            }
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
    {
        addAxes(polygon3);
        f1 = true;
    }

    if (polygon4.vertices.size() != polygon2.vertices.size())
    {
        addAxes(polygon4);
        f2 = true;
    }

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

    cout << m << endl;
    for (int i = 0; i < m; ++i)
    {
        const Vector2D &res = GenSolution(testCases[i]);
        cout << fixed << std::setprecision(5) << res.x << " " << res.y << endl;
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
//         // cout << "Case " << i + 1 << ": " << endl;
//         const Vector2D &res = GenSolution(testCases[i]);
//         cout << fixed << std::setprecision(5) << res.x << " " << res.y << endl;
//         cout.flush();
//     }

//     // send OK after outputting all answer
//     // cout << "OK" << endl;
//     // cout.flush();

//     return 0;
// }