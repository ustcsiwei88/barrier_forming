#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<queue>
#include<unordered_set>
#include<set>
#include<unordered_map>
#include<map>
#include<chrono>
#include<cstdlib>
#include "utils.h"

#include "gurobi_c++.h"

using namespace std;

vector<polygon> obs;
vector<vector<polygon>> s;

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

map<pair<int,int>, int> edge_to_cell;

set<pair<int, int>> straight_lines;

float l, r, u, d;


bool is_collision(const segment & s){
    for(const auto &p: obj){
        for(int i=0;i<p.pts.size();i++){
            int nxt = (i+1) % p.pts.size();
            bool is_eq1 = eq(p.pts[i], s.p1);
            bool is_eq2 = eq(p.pts[i], s.p2);
            bool is_eq3 = eq(p.pts[nxt], s.p1);
            bool is_eq4 = eq(p.pts[nxt], s.p2);
            // if p.pts[i] is between p1 and p2
            if(fabs(cross(p.pts[i] - s.p1, s.p2 - s.p1) ) < 1e-4 && is_between(p.pts[i].x, s.p1.x, s.p2.x) && is_between(p.pts[i].y, s.p1.y, s.p2.y)){
                int pre = (i+p.pts.size()-1) % p.pts.size();
                float theta = is_eq1? (s.p2 - s.p1).atan2(): (s.p1 - s.p2).atan2();
                float theta1 = (p.pts[nxt] - p.pts[i]).atan2();
                float theta2 = (p.pts[pre] - p.pts[i]).atan2();
                if(theta1 > theta2){
                    if(theta > theta1 + 1e-4 || theta < theta2 - 1e-4) {
                        return true;
                    }
                    if(!is_eq1 && !is_eq2){
                        theta = theta + acos(-1.0);
                        if(theta > acos(-1.0)) theta -= 2 * acos(-1.0);
                        if(theta > theta1 + 1e-4 || theta < theta2 - 1e-4) {
                            return true;
                        }
                    }
                }else{
                    if(theta > theta1 + 1e-4 && theta < theta2 - 1e-4) return true;
                    if(!is_eq1 && !is_eq2){
                        theta = theta + acos(-1.0);
                        if(theta > acos(-1.0) * 2) theta -= 2 * acos(-1.0);
                        if(theta > theta1 + 1e-4 && theta < theta2 - 1e-4) return true;
                    }
                }
                
                continue;
            }
            
            double tmp1 = cross(p.pts[i] - s.p1, s.p2 - s.p1) ;
            double tmp2 = cross(p.pts[nxt] - s.p1, s.p2 - s.p1);
            double tmp3 = cross(s.p1 - p.pts[i], p.pts[nxt] - p.pts[i]);
            double tmp4 = cross(s.p2 - p.pts[i], p.pts[nxt] - p.pts[i]);
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
            if(is_between(x, seg.p1.x, seg.p2.x) && is_between(y, seg.p1.y, seg.p2.y) && is_between(x, nseg.p1.x, nseg.p2.x) && is_between(y, nseg.p1.y, nseg.p2.y)){
                pts.emplace_back(x, y);
            }

        }
        if(eq(seg.p1.x, seg.p2.x)) sort(begin(pts), end(pts), [](auto &f, auto &s){return f.y < s.y;});
        else sort(begin(pts), end(pts));
        
        int pre_id = -1;
        for(int i=0; i<pts.size(); i++){
            int id = -1;
            if(i > 0 && eq(pts[i], pts[i-1])) continue;
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
                seg_to_pts[ii].push_back(id);
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

    // int f = 3;
    for(int i = 0; i < points.size(); i++){
        for(int j = 0; j < neighbors[i].size(); j++){
            // cout<<"poly"<<endl;
            vector<int> poly;
            int st = neighbors[i][j];
            int it = st, nit = i;
            // bool convex = true;
            do{
                poly.push_back(it);
                int tmp = nit;
                nit = nxt[nit][it];
                // if(cross(points[tmp] - points[it], points[nit] - points[tmp]) < -1e-3) {break;}
                it = tmp;
                // cout<<"it = "<<it<<endl;
            }while(it != st);
            // if(!convex) continue;
            int min_id = 0;
            for(int k=0; k<poly.size(); k++){
                if(poly[k] < poly[min_id]) min_id = k;
            }
            rotate(begin(poly), begin(poly) + min_id, end(poly));
            if(cell_set.find(poly) == cell_set.end()){
                cell_set.insert(poly);
                int type = -1;
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
                for(int set_id = 0; set_id < s.size(); set_id ++){
                    for(const auto &p: s[set_id]){
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
                            cout << "Set id " << set_id << ' ' << s[set_id].size() << " " << p.pts[0].x << ' ' << p.pts[0].y << endl;
                            
                            type = set_id + 1; break;
                        }
                        
                        // check if the object is inside the cell if it's a point object
                        if(p.pts.size()==1){
                            vec2 v = p.pts[0];
                            int cnt=  0;
                            for(int vn = 0; vn < poly.size(); vn++){
                                vec2 p1 = points[poly[vn]], p2 = points[poly[(vn+1)%poly.size()]];
                                
                                if(p1.y > v.y && p2.y <= v.y || p1.y < v.y && p2.y >= v.y){
                                    if((v.y - p1.y) * (p2.x - p1.x)  / (p2.y - p1.y) + p1.x > v.x) cnt ++;
                                }
                            }
                            if(cnt % 2 == 1){
                                cout << "Set id " << set_id << ' ' << s[set_id].size() << " " << p.pts[0].x << ' ' << p.pts[0].y << endl;
                                type = set_id + 1; break;
                            }else{
                                continue;
                            }
                        }
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
                    type = s.size()+1;
                }
                
                // if(type!=1) cout<< "type = "<<type<<endl;
                if(type == s.size() + 1) continue;

                cell.emplace_back(type, poly);
                cell_edge.emplace_back();
                for(int k =0; k<poly.size(); k++){
                    int fi = poly[k], se = poly[(k+1)%poly.size()];
                    edge_to_cell[make_pair(fi, se)] = cell.size()-1;
                    if(fi > se) swap(fi, se);
                    cell_edge.back().insert(make_pair(fi, se));
                    
                }
            }
        }
    }
    cout<< "Number of cells "<< cell.size()<<endl;
    cout<< "Number of segments candidates " << segments.size()<< ' '<<seg_to_pts.size()<<endl;
    int bits = ceil(log2(s.size()));
    // 0,1 variable for cell
    // 0,1 segment for whole segment
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, true);
    env.start();
    GRBModel model = GRBModel (env);
    vector<vector<GRBVar>> v_c (cell.size(), vector<GRBVar>(bits));
    vector<GRBVar> v_seg (segments.size());
    // cout << "bits = " << bits << endl;
    for(int i=0;i<cell.size();i++){
        for (int j=0; j< bits;j++) v_c[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "v_cell" + to_string(i) + "_" + to_string(j));
    }
    for(int i=0;i<segments.size();i++){
        v_seg[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "v_segment" + to_string(i));
    }
    int cnt = 0;
    GRBLinExpr expr = 0;
    // cout<<"edge_to_cell size "<<edge_to_cell.size()<<endl;
    for(int i=0; i< segments.size(); i++){
        for(int j=1; j<seg_to_pts[i].size(); j++){
            // cout<<"Inner loop"<<endl;
            auto pa = make_pair(seg_to_pts[i][j-1], seg_to_pts[i][j]);
            auto it1 = edge_to_cell.find(pa);
            swap(pa.first, pa.second);
            auto it2 = edge_to_cell.find(pa);
            if(it1 != edge_to_cell.end() && it2 != edge_to_cell.end()){
                int f = it1->second, s = it2->second;
                // if(f==s){cout<<"error"<<endl;}
                for(int j=0; j<bits; j++){
                    model.addConstr(v_seg[i] >= v_c[f][j] - v_c[s][j], "sep_" + to_string(cnt++));
                    model.addConstr(v_seg[i] >= v_c[s][j] - v_c[f][j], "sep_" + to_string(cnt++));
                }
            }
        }
        expr += v_seg[i];
    }
    for(int i=0; i<cell.size() ;i++){
        if(cell[i].first <= s.size() && cell[i].first >=1){
            // cout << "Fixing type " << cell[i].first - 1 << endl;
            for(int j=0; j < bits; j++){
                if((cell[i].first-1) & (1 << j)) model.addConstr(v_c[i][j] >= 1.0, "cell_" + to_string(i));
                else model.addConstr(v_c[i][j] <= 0.0, "cell_" + to_string(i));
            }
        }
    }
    
    model.setObjective(expr, GRB_MINIMIZE);
    model.optimize();
    cout<<"Min number of segments: "<< model.getObjective().getValue()<<endl;
    
    ofstream f_cut("min_segment_cuts.txt");
    for(int i=0;i<segments.size();i++){
        if(v_seg[i].get(GRB_DoubleAttr_X) == 1.0){
            // cout<< "whywhy"<<endl;
            f_cut << segments[i].p1.x<<' ' << segments[i].p1.y << ' ' << segments[i].p2.x <<' ' << segments[i].p2.y<<endl;
        }
    }
    f_cut.close();

    ofstream f_cells("cells.txt");
    for(int cnt = 0; cnt < cell.size(); cnt ++){
        auto tmp = cell[cnt].second;
        f_cells << cell[cnt].first << ' ';
        for(int i:tmp){
            // f_cells <
            f_cells << points[i].x << ' '<< points[i].y << ' ' ;
        }
        f_cells << endl;
    }
    f_cells.close();
}

void extend(segment &seg){
    vector<vec2> ext_pts{seg.p1, seg.p2};
    polygon tmp;
    tmp.pts = vector<vec2>{vec2(l, d), vec2(r, d), vec2(r, u), vec2(l, u)};
    obj.push_back(tmp);
    for (const auto o: obj){
        for(int ei = 0; ei< o.pts.size(); ei++){
            int ni = (ei + 1)%o.pts.size();
            vec2 tmp;
            if(get_intersection(seg, segment(o.pts[ei], o.pts[ni]), tmp)){
                if(!is_between(tmp.x, o.pts[ei].x, o.pts[ni].x) || !is_between(tmp.y, o.pts[ei].y, o.pts[ni].y)) continue;
                bool iseq = false;
                for(const auto& pt: ext_pts){
                    if(eq(pt, tmp)){
                        iseq = true; break;
                    }
                }
                if(!iseq) ext_pts.push_back(tmp);
            }
        }
    }
    obj.pop_back();
    sort(begin(ext_pts), end(ext_pts));
    if(seg.p2 < seg.p1) swap(seg.p1, seg.p2);
    auto ptr = lower_bound(begin(ext_pts), end(ext_pts), seg.p1);
    while(ptr!=ext_pts.begin() && !is_collision(segment(seg.p2, *(ptr-1)))) seg.p1 = *(--ptr);
    ptr = lower_bound(begin(ext_pts), end(ext_pts), seg.p2);
    while(ptr != ext_pts.end() && !is_collision(segment(seg.p1, *(ptr)))) seg.p2 = *ptr ++;
}

int main(){
    ifstream f("input.txt");
    ofstream f_seg("segments.txt");
    int obs_num, n_pts, s_num;
    float x, y;
    f >> l >> d >> r >> u >> s_num;
    // vector<int> set_nums;
    segments.emplace_back(vec2(l, d), vec2(r, d));
    segments.emplace_back(vec2(r, d), vec2(r, u));
    segments.emplace_back(vec2(r, u), vec2(l, u));
    segments.emplace_back(vec2(l, u), vec2(l, d));

    polygon workspace;
    workspace.pts = vector<vec2>{vec2(l, d), vec2(r, d), vec2(r, u), vec2(l, u)};
    // obj.push_back(workspace);
    vector<vector<polygon>*> input{&obs};
    s.resize(s_num);
    for(int i=0; i < s_num; i++){
        input.push_back(&s[i]);
    }

    for(auto p: input){
        vector<polygon>& o = *p;
        f >> obs_num;
        // cout<<"obs_num" << obs_num<<endl;
        for(int i=0; i<obs_num; i++){
            o.emplace_back();
            f >> n_pts;

            for(int j=0; j<n_pts; j++){
                f >> x >> y;
                cout<<x<<' '<<y<<endl;
                o.back().pts.emplace_back(x, y);
            }
            obj.push_back(o.back());
        }
    }

using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    for(int i=0;i<obj.size();i++){
        for(int j=0;j<obj[i].pts.size();j++){
            // point to point segment
            vec2 p = obj[i].pts[j];
            for(int k=j+1;k<obj[i].pts.size();k++){
                segment seg(p, obj[i].pts[k]);
                if(!is_collision(seg)){
                    // extend to both sides
                    extend(seg);
                    segments.push_back(seg);
                }
            }
            for(int k = i+1; k< obj.size(); k++){
                for(int l=0;l<obj[k].pts.size();l++){
                    auto q = obj[k].pts[l];
                    auto v = p - obj[k].pts[l];
                    swap(v.x, v.y);
                    v.x = -v.x;
                    double val = v.norm();
                    v.x /= val * 1e2; v.y /= val * 1e2;
                    if(obj[i].pts.size() == 1){
                        if(obj[k].pts.size()== 1){
                            for(auto seg: vector<segment>{segment(p + v, q + v), segment(p + v, q - v), segment(p - v, q + v),
                                    segment(p - v, q - v)}){
                                if(!is_collision(seg)){
                                    extend(seg);
                                    segments.push_back(seg);
                                }
                            }
                        }else{
                            for(auto seg: vector<segment>{segment(p + v, q), segment(p - v, q)}){
                                if(!is_collision(seg)){
                                    extend(seg);
                                    segments.push_back(seg);
                                }
                            }
                        }
                    }else{
                        if(obj[k].pts.size() == 1){
                            for(auto seg: vector<segment>{segment(p, q - v), segment(p, q + v)}){
                                if(!is_collision(seg)){
                                    extend(seg);
                                    segments.push_back(seg);
                                }
                            }
                        }else{
                            segment seg(p, obj[k].pts[l]);
                            if(!is_collision(seg)){
                                extend(seg);
                                segments.push_back(seg);
                            }
                        }
                    }
                }
            }
            
        }
    }
    // remove redundancy
    for(auto & seg: segments){
        if(seg.p2 < seg.p1) swap(seg.p1, seg.p2);
    }
    sort(begin(segments), end(segments));
    for(int i=0, j=1; j < segments.size();){
        bool equals = false;
        for (int k=0; k<=i; k++){
            if(eq(segments[k].p1, segments[j].p1) && eq(segments[k].p2, segments[j].p2)) {
                equals = true;
                break;
            }
        }
        if(equals) j++;
        else segments[++i] = segments[j++];
        if(j == segments.size()){
            segments.resize(i+1); break;
        }
    }
    for(const auto &seg: segments){
        f_seg << seg.p1.x <<' '<< seg.p1.y << ' '<<seg.p2.x <<' '<<seg.p2.y<<endl;
    }
    cell_construction();

    auto timecost = duration_cast<duration<double>>(high_resolution_clock::now() - start_time).count();

    fstream file;
    file.open("time.txt", std::ios_base::app);
    file << timecost << endl;
    file.close();
    return 0;
}

