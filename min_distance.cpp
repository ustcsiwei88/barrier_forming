#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<queue>
#include<unordered_set>
#include<set>
#include<unordered_map>
#include "utils.h"

using namespace std;

vector<polygon> obs;
vector<polygon> s1;
vector<polygon> s2;

vector<polygon> obj;

vector<segment> segments;

vector<vector<int>> seg_to_pts;
vector<segment> edges;

vector<vec2> points;
vector<vector<int>> neighbors;
vector<unordered_map<int,int>> nxt;
// unordered_set<vec2> set_of_points;

// 0 obstacle 1 free space 2 start 3 goal
vector<pair<int, vector<int> > > cell;
set<vector<int>> cell_set;
vector<set<pair<int,int>>> cell_edge; 

float l, r, u, d;


namespace maxflow{
    int num;
    int s, t;
    vector<vector<float>> capa;
    vector<vector<int>> nei;
    // vector<vector<bool>> b_nei;
    void init(){
        capa.resize(num, vector<float>(num, 0.0));
        nei.resize(num);
        cout<<"Init succ"<<endl;
    }
    void add_edge(int s, int t, float dis){
        nei[s].push_back(t); nei[t].push_back(s);
        capa[s][t] = dis; capa[t][s] = dis;
    }
    vector<pair<int,int>> cut;
    float compute_maxflow(){
        float res = 0;
        vector<int> parent(num, -1);
        while(1){
            fill(begin(parent), end(parent), -1);
            queue<int> q;
            q.push(s);
            parent[s] = s;
            bool updated = false;
            while(q.size()){
                // cout<<"Updating"<<endl;
                int i = q.front(); q.pop();
                for(int j: nei[i]){
                    if(parent[j] ==-1 && capa[i][j] > 1e-6){
                        parent[j] = i;
                        q.push(j);
                        if(j == t){
                            updated = true;
                            float gain = 1e9;
                                cout<<i<<' '<<j<< ' '<<s<<' '<<t<<endl;
                            for(j = t, i = parent[t]; j!=s; j = i, i=parent[i]){
                                gain = min(capa[i][j], gain);
                                // cout<<i<<' '<<j<<endl;
                            }
                                // cout<<i<<' '<<j<< ' '<<s<<' '<<t<<endl;
                            res += gain;
                            cout<<"Gained "<<gain<<endl;
                            for(j = t, i = parent[t]; j!=s; j = i, i=parent[i]){
                                capa[i][j] -= gain;
                                capa[j][i] += gain;
                            }
                            updated = true; goto r_end;
                        }
                    }
                }
            }
            r_end:
            if(!updated){
                break;
            }
        }
        cut.clear();
        for(int i=0;i<parent.size();i++){
            if(parent[i] != -1){
                for(int j:nei[i]){
                    if(parent[j] == -1){
                        cut.emplace_back(i, j);
                    }
                }
            }
        }
        return res;
    }
};


// Compute line segments and do plane tessellation

// no need
// bool is_inside(const vec2 &pt, const polygon &p){
//     int cnt = 0;
//     for(int i=0; i<p.pts.size();i++){
//         vec2 p1 = p.pts[i];
//         int nxt = (i+1) % p.pts.size();
//         vec2 p2 = p.pts[nxt];
        
//     }
// }



bool is_collision(const segment & s){
    // cout<<"s="<<s.p1.x <<' '<<s.p1.y<< ' '<<s.p2.x<<' '<<s.p2.y<<endl;
    for(const auto &p: obj){
        for(int i=0;i<p.pts.size();i++){
            // cout<<i<<endl;
            int nxt = (i+1) % p.pts.size();
            bool is_eq1 = eq(p.pts[i], s.p1);
            bool is_eq2 = eq(p.pts[i], s.p2);
            bool is_eq3 = eq(p.pts[nxt], s.p1);
            bool is_eq4 = eq(p.pts[nxt], s.p2);
            if(fabs(cross(p.pts[i] - s.p1, s.p2 - s.p1) ) < 1e-6 && is_between(p.pts[i].x, s.p1.x, s.p2.x) && is_between(p.pts[i].y, s.p1.y, s.p2.y)){
                int pre = (i+p.pts.size()-1) % p.pts.size();
                float theta = is_eq1? (s.p2 - s.p1).atan2(): (s.p1 - s.p2).atan2();
                float theta1 = (p.pts[nxt] - p.pts[i]).atan2();
                float theta2 = (p.pts[pre] - p.pts[i]).atan2();
                if(theta1 > theta2){
                    if(theta > theta1 + 1e-6 || theta < theta2 - 1e-6) {
                        // cout<<theta<<' '<<theta1 <<' '<<theta2<<"==================="<<endl;
                        return true;
                    }
                    if(!is_eq1 && !is_eq2){
                        theta = theta + acos(-1.0);
                        if(theta > acos(-1.0)) theta -= 2 * acos(-1.0);
                        if(theta > theta1 + 1e-6 || theta < theta2 - 1e-6) {
                            return true;
                        }
                    }
                }else{
                    if(theta > theta1 + 1e-6 && theta < theta2 - 1e-6) return true;
                    if(!is_eq1 && !is_eq2){
                        theta = theta + acos(-1.0);
                        if(theta > acos(-1.0) * 2) theta -= 2 * acos(-1.0);
                        if(theta > theta1 + 1e-6 && theta < theta2 - 1e-6) return true;
                    }
                }
                
                continue;
            }
            
            float tmp1 = cross(p.pts[i] - s.p1, s.p2 - s.p1) ;
            float tmp2 = cross(p.pts[nxt] - s.p1, s.p2 - s.p1);
            float tmp3 = cross(s.p1 - p.pts[i], p.pts[nxt] - p.pts[i]);
            float tmp4 = cross(s.p2 - p.pts[i], p.pts[nxt] - p.pts[i]);
            if(tmp1 * tmp2 < -1e-6 && tmp3 * tmp4 < -1e-6) return true;
            if(tmp1 * tmp2 > 1e-6 || tmp3 * tmp4 > 1e-6) continue;
            if(tmp3 * tmp4 > -1e-6 && !is_eq1 && !is_eq2 && !is_eq3 && !is_eq4){
                if(tmp3 < -1e-6 || tmp4 < -1e-6) return true;
            }
        }
    }
    return false;
}



void cell_construction(){
    seg_to_pts.resize(segments.size());
    for(int ii=0; ii<segments.size();ii++){
        segment seg = segments[ii];
        vector<vec2> pts{seg.p1, seg.p2};
        for(int j=0; j<segments.size();j++){
            if(ii==j) continue;
            segment nseg = segments[j];
            float a1 = seg.p2.y - seg.p1.y;
            float b1 = seg.p1.x - seg.p2.x;
            float c1= cross(seg.p2, seg.p1);
            float a2 = nseg.p2.y - nseg.p1.y;
            float b2 = nseg.p1.x - nseg.p2.x;
            float c2= cross(nseg.p2, nseg.p1);
            float deter = a1 * b2 - a2 * b1;
            if(deter > -1e-6 && deter < 1e-6) continue;
            float x = (b1 * c2 - b2 * c1) / deter;
            float y = (a2 * c1 - a1 * c2) / deter;
            // if(fabs(y) <= 1e-6 && x > 0){
            //     cout<<"Got "<<x<<' '<<y<<endl;
            //     cout<<seg.p1<< seg.p2 << nseg.p1<<nseg.p2<<endl;
            // }
            if(is_between(x, seg.p1.x, seg.p2.x) && is_between(y, seg.p1.y, seg.p2.y) && is_between(x, nseg.p1.x, nseg.p2.x) && is_between(y, nseg.p1.y, nseg.p2.y)){
                pts.emplace_back(x, y);
                // if(fabs(y) <= 1e-6 && x > 0)cout<<"Got "<<x<<' '<<y<<endl;
            }

        }
        sort(begin(pts), end(pts));
        // cout<<"List: ";
        // for(const auto& pt: pts) cout<<pt<<' '; cout<<endl;
        // cout<<"intersection points constructed"<<endl;
        int pre_id = -1;
        for(int i=0; i<pts.size(); i++){
            int id = -1;
            for(int j = 0; j< points.size();j++){
                if(eq(points[j], pts[i])){
                    id = j; break;
                }
            }
            if(id == -1) {
                id = points.size();
                points.push_back(pts[i]);
                neighbors.emplace_back();
            }
            if(pre_id == -1 || pre_id == id);
            else{
                neighbors[pre_id].push_back(id);
                neighbors[id].push_back(pre_id);
            }
            if(pre_id != id){
                seg_to_pts[i].push_back(id);
            }
            pre_id = id;
        }
    }
    // cout<<"HERE"<<endl;
    // cout<<points.size() << ' '<<neighbors.size()<<endl;
    nxt.resize(neighbors.size());
    for(int i =0; i< points.size(); i++){
        sort(begin(neighbors[i]), end(neighbors[i]));
        neighbors[i].resize(unique(begin(neighbors[i]), end(neighbors[i]) ) - begin(neighbors[i]));
        // vec2 p = points[i];
        // cout<<"HERE"<<endl;
        sort(begin(neighbors[i]), end(neighbors[i]), [&](const auto &j, const auto &k){
            if(fabs((points[j] - points[i]).atan2() - (points[k] - points[i]).atan2()) < 1e-5){
                cout<<"error occurred"<<endl;
                cout<<points[j]<<' '<<points[i]<<' '<<points[k]<<endl;
            }
            return (points[j] - points[i]).atan2() < (points[k] - points[i]).atan2();
        });
        // cout<<"HERE"<<endl;
        for(int j=0; j<neighbors[i].size(); j++){
            nxt[i][neighbors[i][j]] = neighbors[i][(j + neighbors[i].size() - 1) % neighbors[i].size()];
        }
    }

    int f = 3;
    for(int i = 0; i < points.size(); i++){
        for(int j = 0; j < neighbors[i].size(); j++){
            // cout<<"poly"<<endl;
            vector<int> poly;
            int st = neighbors[i][j];
            int it = st, nit = i;
            do{
                poly.push_back(it);
                int tmp = nit;
                nit = nxt[nit][it];
                it = tmp;
                // cout<<"it = "<<it<<endl;
            }while(it != st);
            int min_id = 0;
            for(int k=0; k<poly.size(); k++){
                if(poly[k] < poly[min_id]) min_id = k;
            }
            rotate(begin(poly), begin(poly) + min_id, end(poly));
            if(cell_set.find(poly) == cell_set.end()){
                cell_set.insert(poly);
                int type = 1;
                for(const auto &p: obs){
                    bool is_in = true;
                    int cnt = 0;
                    // cout<<"HERE"<<endl;
                    for(const auto &v: p.pts){
                        bool is_in2 = false;
                        for(int vn: poly){
                            if(eq(points[vn], v)){
                                is_in2 = true; break;
                            }
                        }
                        if(!is_in2){
                            is_in = false;
                            break;
                        }cnt++;
                    }
                    // cout<<"matched "<<cnt<<"points"<<endl;
                    if(is_in){
                        type = 0; break;
                    }
                }
                for(const auto &p: s1){
                    bool is_in = true;
                    for(const auto &v: p.pts){
                        bool is_in2 = false;
                        for(int vn: poly){
                            if(eq(points[vn], v)){
                                is_in2 = true; break;
                            }
                        }
                        if(!is_in2){
                            is_in = false;
                            break;
                        }
                    }
                    if(is_in){
                        type = 2; break;
                    }
                }
                for(const auto &p: s2){
                    bool is_in = true;
                    for(const auto &v: p.pts){
                        bool is_in2 = false;
                        for(int vn: poly){
                            if(eq(points[vn], v)){
                                is_in2 = true; break;
                            }
                        }
                        if(!is_in2){
                            is_in = false;
                            break;
                        }
                    }
                    if(is_in){
                        type = 3; break;
                    }
                }
                bool is_in = true;
                for(const auto &v: vector<vec2>{vec2(l,d), vec2(r,d), vec2(r, u), vec2(l,u)}){
                    bool is_in2 = false;
                    for(int vn: poly){
                        if(eq(points[vn], v)){
                            is_in2 = true; break;
                        }
                    }
                    if(!is_in2){
                        is_in = false;
                        break;
                    }
                }
                if(is_in){
                    type = 4;
                }
                
                // if(type!=1) cout<< "type = "<<type<<endl;
                if(type == 0 || type == 4) continue;

                cell.emplace_back(type, poly);
                cell_edge.emplace_back();
                for(int k =0; k<poly.size(); k++){
                    int fi = poly[k], se = poly[(k+1)%poly.size()];
                    if(fi > se) swap(fi, se);
                    cell_edge.back().insert(make_pair(fi, se));
                }
                // if(f == 0){
                //     for(int vn: poly){
                //         cout<< points[vn];
                //     }cout<<endl;
                    
                // }f--;
            }
        }
    }
    cout<< "Number of cells "<< cell.size()<<endl;
    using namespace maxflow;
    num = cell.size() + 2;


    // min_distance

    init();
    s = num - 2, t = num - 1;
    // cout<<"edge size "<< (cell_edge.size())<<endl;
    for(int i=0;i<cell.size();i++){
        if(cell[i].first == 2){
            add_edge(s, i, 1e9);
        }
        if(cell[i].first == 3){
            add_edge(t, i, 1e9);
        }
        for(int j=i+1;j<cell.size();j++){
            float total_dis = .0;
            for(const auto & e: cell_edge[i]){
                if(cell_edge[j].find(e) != cell_edge[j].end()){
                    // if(e.first >= points.size() || e.second >= points.size()){
                    //     // cout<<"Error here"<< e.first<<' '<<e.second<<endl;
                    // }
                    total_dis += (points[e.first] - points[e.second]).norm();
                }
            }
            if(total_dis > 1e-6) {
                add_edge(i, j, total_dis);
            }
        }
        // cout<<"FFF "<<i<<endl;
    }
    cout<<"min distance " << compute_maxflow()<<endl;
    ofstream f_cut("cuts.txt");
    for(const auto [i, j]: cut){
        for(const auto & e: cell_edge[i]){
            if(cell_edge[j].find(e) != cell_edge[j].end()){
                f_cut<< points[e.first].x <<' ' << points[e.first].y <<' '<<points[e.second].x <<' ' << points[e.second].y <<endl;
            }
        }
    }
}

int main(){
    ifstream f("input.txt");
    ofstream f_seg("segments.txt");
    int obs_num, s1_num, s2_num, n_pts;
    float x, y;
    f >> l >> d >> r >> u;
    segments.emplace_back(vec2(l, d), vec2(r, d));
    segments.emplace_back(vec2(r, d), vec2(r, u));
    segments.emplace_back(vec2(r, u), vec2(l, u));
    segments.emplace_back(vec2(l, u), vec2(l, d));
    for(auto p: vector<vector<polygon>*>{&obs, &s1, &s2}){
        vector<polygon>& o = *p;
        f >> obs_num;
        // cout<<"obs_num" << obs_num<<endl;
        for(int i=0; i<obs_num; i++){
            o.emplace_back();
            f >> n_pts;
            cout<<n_pts<<endl;
            for(int j=0; j<n_pts; j++){
                f >> x >> y;
                cout<<x<<' '<<y<<endl;
                o.back().pts.emplace_back(x, y);
            }
            obj.push_back(o.back());
        }
    }
    for(int i=0;i<obj.size();i++){
        for(int j=0;j<obj[i].pts.size();j++){
            // point to point segment
            vec2 p = obj[i].pts[j];
            for(int k=j+1;k<obj[i].pts.size();k++){
                segment seg(p, obj[i].pts[k]);
                if(!is_collision(seg)){
                    segments.push_back(seg);
                }
            }
            for(int k = i+1; k< obj.size(); k++){
                for(int l=0;l<obj[k].pts.size();l++){
                    segment seg(p, obj[k].pts[l]);
                    if(!is_collision(seg)){
                        segments.push_back(seg);
                    }
                }
            }
            // point to line segment
            for(int k = 0; k< obj.size(); k++){
                if(k == i) continue;
                for(int l=0; l<obj[k].pts.size(); l++){
                    vec2 p1 = obj[k].pts[l];
                    int nxt = (l+1) % obj[k].pts.size();
                    vec2 p2 = obj[k].pts[nxt];
                    float A = p2.y - p1.y, B = p1.x - p2.x, C = cross(p2, p1);
                    float lambda = (A * p.x + B * p.y + C ) / (A*A + B*B);
                    vec2 pedal(p.x - lambda * A, p.y - lambda * B);
                    if(pedal.x > min(p1.x, p2.x) && pedal.x <max(p1.x, p2.x) && pedal.y > min(p1.y, p2.y) && pedal.y < max(p1.y, p2.y)){
                        segment seg(p, pedal);
                        if(!is_collision(seg)){
                            segments.push_back(seg);
                        }
                    }
                }
            }
            // point to boundary segment
            for(const auto &seg: vector<segment>{segment(p, vec2(l, p.y)), segment(p, vec2(r, p.y)), segment(p, vec2(p.x, d)), segment(p, vec2(p.x, u))}){
                if(!is_collision(seg)){
                    segments.push_back(seg);
                }
            }
        }
    }
    // remove redundancy?
    // for(int i=0; i<segments.size();i++){
    //     segment seg = segments[i];
    //     if(seg.p1.x > seg.p2.x){
    //         swap(seg.p1, seg.p2);
    //     }else if(seg.p1.x == seg.p2.x && seg.p1.y > seg.p2.y){
    //         swap(seg.p1, seg.p2);
    //     }
    //     for(int j=0; j<segments.size(); j++){
    //         if(i==j) continue;
    //         segment nseg = segments[j];
    //         if()
    //     }
    // }
    for(const auto &seg: segments){
        f_seg << seg.p1.x <<' '<< seg.p1.y << ' '<<seg.p2.x <<' '<<seg.p2.y<<endl;
    }
    cell_construction();
    return 0;
}