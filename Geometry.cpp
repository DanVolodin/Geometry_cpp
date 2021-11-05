#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

typedef long double ld;
typedef long long ll;

const ld pi = acos(-1), eps = 1e-6;
const int INF = 1e9 + 7;


struct Point {
    ld x, y;
    Point() {}
    Point(ld x, ld y) : x(x), y(y) {}

    bool operator ==(Point b) {
        return x == b.x && y == b.y;
    }

    friend istream& operator >>(istream& in, Point& a) {
        in >> a.x >> a.y;
        return in;
    }

    friend ostream& operator <<(ostream& out, Point& a) {
        out << a.x << " " << a.y;
        return out;
    }
};

struct Vector {
    ld x, y;
    Vector() {}
    Vector(ld x1, ld y1, ld x2, ld y2) : x(x2 - x1), y(y2 - y1) {}
    Vector(ld x, ld y) : x(x), y(y) {}
    Vector(Point a, Point b) : x(b.x - a.x), y(b.y - a.y) {}

    Vector operator +(Vector& b) {
        return Vector(x + b.x, y + b.y);
    }

    Vector operator -(Vector& b) {
        return Vector(x - b.x, y - b.y);
    }

    ld operator *(Vector& b) {
        return x * b.x + y * b.y;
    }

    ld operator /(Vector b) {
        return x * b.y - y * b.x;
    }

    friend istream& operator >>(istream& in, Vector& a) {
        in >> a.x >> a.y;
        return in;
    }

    friend ostream& operator <<(ostream& out, Vector& a) {
        out << a.x << " " << a.y;
        return out;
    }
};

struct Polygon {
    int sz;
    vector <Point> vertex;
    Polygon() {}
    Polygon(int sz_) {
        sz = sz_;
        vertex.resize(sz);
    }
    Polygon(int sz_, vector <Point> points) {
        sz = sz_;
        vertex.resize(sz);
        copy(points.begin(), points.end(), vertex.begin());
    }

    friend istream& operator >>(istream& in, Polygon& p) {
        in >> p.sz;
        p.vertex.resize(p.sz);
        for (auto& A : p.vertex) in >> A;
        return in;
    }

    friend ostream& operator <<(ostream& out, Polygon& p) {
        out << p.sz << "\n";
        for (auto& A : p.vertex) out << A << "\n";
        return out;
    }
};

struct Line {
    ld A, B, C;
    Line() {};
    Line(ld A, ld B, ld C) : A(A), B(B), C(C) {}
    Line(Point A1, Point A2) {
        A = A2.y - A1.y;
        B = A1.x - A2.x;
        C = A2.x * A1.y - A1.x * A2.y;
    }
    Line(Point P, Vector V, int type) {
        if (type == 0) { //parallel
            A = V.y;
            B = -V.x;
            C = V.x * P.y - V.y * P.x;
        }
        else if (type == 1){ //perpendicular
            A = V.x;
            B = V.y;
            C = -(V.x * P.x + V.y * P.y);
        }
    }

    friend istream& operator >>(istream& in, Line& a) {
        in >> a.A >> a.B >> a.C;
        return in;
    }

    friend ostream& operator <<(ostream& out, Line& a) {
        out << a.A << " " << a.B << " " << a.C;
        return out;
    }
};

Point tmp;

bool point_on_ray(Point P, Point O, Point A) {
    Vector op(O, P), oa(O, A);
    return op / oa == 0 && op * oa >= 0;
}

bool point_on_segment(Point P, Point A, Point B) {
    return point_on_ray(P, A, B) && point_on_ray(P, B, A);
}

double len(Vector a) {
    return hypot(a.x, a.y);
}

double dist_point_segment(Point P, Point A, Point B) {
    if (point_on_segment(P, A, B)) return 0;
    Vector ap(A, P), ab(A, B);
    Vector bp(B, P), ba(B, A);
    if (ap * ab < 0) return len(ap);
    if (bp * ba < 0) return len(bp);
    return abs(ap / ab) / len(ab);
}

double angle(Vector a, Vector b) {
    return abs(atan2(a / b, a * b));
}

double angle(Point A, Point O, Point B) {
    return angle(Vector(O, A), Vector(O, B));
}

bool diff(ld x, ld y) {
    return x <= 0 && y >= 0 || x >= 0 && y <= 0;
}

bool point_in_angle(Point P, Point O, Vector a, Vector b) {
    if (P == O) return true;
    Vector p(O, P);
    double aop = angle(a, p);
    double bop = angle(b, p);
    double aob = angle(a, b);
    if (abs(aop + bop - aob) > eps) return false;
    if (!diff(p / a, p / b)) return false;
    return true;
}

bool intersect(ld a, ld b, ld c, ld d) {
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    return min(b, d) >= max(a, c);
}

bool segments_intersect(Point a, Point b, Point c, Point d) {
    Vector ab(a, b), ba(b, a), cd(c, d), dc(d, c);
    Vector ca(c, a), cb(c, b), da(d, a), db(d, b);
    Vector ac(a, c), ad(a, d), bc(b, c), bd(b, d);
    return diff(cd / ca, cd / cb) && diff(dc / da, dc / db) &&
    diff(ab / ac, ab / ad) && diff(ba / bc, ba / bd) &&
    intersect(a.x, b.x, c.x, d.x) && intersect(a.y, b.y, c.y, d.y);
}

Point intersect_p(Line a, Line b) {
    return Point(-(a.B * b.C - b.B * a.C) * 1.0 / (a.B * b.A - b.B * a.A),
                 -(a.A * b.C - b.A * a.C) * 1.0 / (a.A * b.B - b.A * a.B));
}

bool is_convex(Polygon p) {
    Point a = p.vertex[p.sz - 1];
    Point b = p.vertex[0];
    Point c = p.vertex[1];
    ld base = Vector(a, b) / Vector(b, c);
    for (int i = 0; i < p.sz - 2; i++) {
        a = p.vertex[i];
        b = p.vertex[i + 1];
        c = p.vertex[i + 2];
        ld cur = Vector(a, b) / Vector(b, c);
        if (diff(base, cur)) return false;
    }
    return true;
}

bool point_in_polygon(Point A, Polygon p) {
    double sum = 0;
    for (int i = 0; i < p.sz; i++) {
        if (point_on_segment(A, p.vertex[i], p.vertex[(i + 1) % p.sz])) return true;
        Vector a(A, p.vertex[i]), b(A, p.vertex[(i + 1) % p.sz]);
        sum += angle(a, b);
    }
    return abs(abs(sum) - 2 * pi) < eps;
}

double area(Polygon p) {
    ld res = 0;
    for (int i = 0; i < p.sz; i++) {
        Point A = p.vertex[i], B = p.vertex[(i + 1) % p.sz];
        res += (A.x * B.y - A.y * B.x);
    }
    return abs(res);
}

Vector normalize(Vector a, double k) {
    ld d = hypot(a.x, a.y);
    return Vector(k * a.x / d, k * a.y / d);
}

Line bisector(Point A, Point O, Point B) {
    if (abs(angle(A, O, B) - pi) < eps) {
        Line a(O, A);
        return Line(O, Vector(a.A, a.B), 0);
    }
    Vector oa(O, A), ob(O, B);
    Vector oa_ = normalize(oa, len(ob));
    Point A_(O.x + oa_.x, O.y + oa_.y);
    Line a(O, A_), b(O, B);
    Vector na(a.A, a.B), nb(b.A, b.B);
    Line b_perp(B, nb, 0), a_perp(A_, na, 0);
    Point C = intersect_p(a_perp, b_perp);
    return Line(O, C);
}

bool cmp(Point& a, Point& b) {
    Vector oa(tmp, a), ob(tmp, b);
    if (oa / ob == 0) return len(oa) < len(ob);
    return oa / ob > 0;
}

Polygon convex_hull(vector <Point>& points) {
    Point st(INF, INF);
    for (auto A : points) {
        if (A.x < st.x || A.x == st.x && A.y < st.y) st = A;
    }
    tmp = st;
    sort(points.begin(), points.end(), cmp);
    vector <Point> ans;
    ans.push_back(st);
    int sz = 1;
    for (int i = 1; i < points.size(); i++) {
        while (sz > 0 && ans[sz - 1] == points[i] ||
        sz > 1 && Vector(ans[sz - 2], ans[sz - 1]) / Vector(ans[sz - 1], points[i]) <= 0) {
            ans.pop_back();
            sz--;
        }
        sz++;
        ans.push_back(points[i]);
    }
    return Polygon(sz, ans);
}


int main(){
    /*To solve geometry problems write your code here using functions above*/
    return 1;
}
