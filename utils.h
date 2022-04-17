#include<vector>
#include<iostream>
#include<algorithm>
#include<cmath>

#define eps 1e-6

using namespace std;
struct vec2{
    double x, y;
    vec2(float x=0.0, float y=0.0):x(x), y(y){}
    double atan2(){return std::atan2(y, x);}
    double norm(){return hypot(x, y);}
};

vec2 operator + (const vec2& v1, const vec2& v2){
    return vec2(v1.x + v2.x, v1.y + v2.y);
}

vec2 operator - (const vec2& v1, const vec2& v2){
    return vec2(v1.x - v2.x, v1.y - v2.y);
}

bool operator < (const vec2& v1, const vec2 &v2){
    if(fabs(v1.x - v2.x) > 1e-6){
        if(v1.x < v2.x) return true;
        if(v1.x > v2.x) return false;
    }if(v1.y < v2.y) return true;
    return false;
}



ostream& operator <<(ostream &os, const vec2 & v){
    os << v.x <<','<<v.y<<' ';
    return os;
}

float cross(const vec2& v1, const vec2& v2){
    return v1.x * v2.y - v1.y * v2.x;
}

struct polygon{
    vector<vec2> pts;
};

struct segment{
    vec2 p1, p2;
    segment(){}
    segment(const vec2 &p1, const vec2 &p2):p1(p1), p2(p2){}
};

bool operator < (const segment& s1, const segment & s2){
    if(s1.p1 < s2.p1) return true;
    if(s2.p1 < s1.p1) return false;
    if(s1.p2 < s2.p2) return true;
    return false;
}

inline bool is_between(float a, float b, float c){
    return a <= max(b,c) + 1e-6 && a >= min(b,c) - 1e-6;
}

inline bool eq(const vec2 &p1, const vec2 &p2){
    return fabs(p1.x - p2.x) + fabs(p1.y - p2.y) < 1e-4;
}

inline bool eq(const float a, const float b){
    return fabs(a - b) < 1e-5;
}

bool operator == (const segment& s1, const segment & s2){
    return eq(s1.p1, s2.p1) && eq(s1.p2, s2.p2);
}

bool get_intersection(const segment& seg, const segment& nseg, vec2& res){
    float a1 = seg.p2.y - seg.p1.y;
    float b1 = seg.p1.x - seg.p2.x;
    float c1= cross(seg.p2, seg.p1);
    float a2 = nseg.p2.y - nseg.p1.y;
    float b2 = nseg.p1.x - nseg.p2.x;
    float c2= cross(nseg.p2, nseg.p1);
    float deter = a1 * b2 - a2 * b1;
    if(deter > -1e-6 && deter < 1e-6) return false;
    res.x = (b1 * c2 - b2 * c1) / deter;
    res.y = (a2 * c1 - a1 * c2) / deter;
    return true;
}
