#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<queue>
#include<unordered_set>
#include<set>
#include<unordered_map>
#include<map>
#include "utils.h"

#include "gurobi_c++.h"

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

map<pair<int,int>, int> edge_to_cell;

set<pair<int, int>> straight_lines;

float l, r, u, d;


namespace maxflow{
    int num;
    int s, t;
    vector<vector<float>> capa;
    vector<vector<int>> nei;
    // vector<vector<bool>> b_nei;
    vector<bool> visited;
    void init(){
        capa.resize(num, vector<float>(num, 0.0));
        nei.resize(num);
        visited.resize(num);
    }
    void add_edge(int s, int t, float dis){
        nei[s].push_back(t); nei[t].push_back(s);
        capa[s][t] = dis; //capa[t][s] = dis;
    }
    // vector<pair<int,int>> cut;
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
                            for(j = t, i = parent[t]; j!=s; j = i, i=parent[i]){
                                gain = min(capa[i][j], gain);
                            }
                            res += gain;
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
        fill(begin(visited), end(visited), false);
        for(int i=0;i<parent.size();i++){
            if(parent[i] != -1){
                visited[i] = true;
            }
        }

        return res;
    }
};




bool is_collision(const segment & s){
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
            if(is_between(x, seg.p1.x, seg.p2.x) && is_between(y, seg.p1.y, seg.p2.y) && is_between(x, nseg.p1.x, nseg.p2.x) && is_between(y, nseg.p1.y, nseg.p2.y)){
                pts.emplace_back(x, y);
            }

        }
        sort(begin(pts), end(pts));
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
                seg_to_pts[ii].push_back(id);
            }
            pre_id = id;
        }
    }
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
        for(int j=0; j<neighbors[i].size(); j++){
            nxt[i][neighbors[i][j]] = neighbors[i][(j + neighbors[i].size() - 1) % neighbors[i].size()];
        }
    }

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
                if(/*type == 0 ||*/ type == 4) continue;

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
    
    // 0,1 variable for cell
    // 0,1 variable for unit segment
    // 0,1 segment for whole segment
    GRBEnv env = GRBEnv(true);
    env.set(GRB_IntParam_OutputFlag, true);
    env.start();
    GRBModel model = GRBModel (env);
    GRBVar *v_c = new GRBVar[cell.size()];
    GRBVar *v_seg = new GRBVar[segments.size()];
    for(int i=0;i<cell.size();i++){
        v_c[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "v_cell" + to_string(i));
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
                if(f==s){cout<<"error"<<endl;}
                model.addConstr(v_seg[i] >= v_c[f] - v_c[s], "sep_" + to_string(cnt++));
                model.addConstr(v_seg[i] >= v_c[s] - v_c[f], "sep_" + to_string(cnt++));
            }
        }
        expr += v_seg[i];
    }
    for(int i=0; i<cell.size() ;i++){
        if(cell[i].first == 2){
            model.addConstr(v_c[i] >= 1.0, "cell_" + to_string(i));
        }else if (cell[i].first == 3){
            model.addConstr(v_c[i] <= 0.0, "cell_" + to_string(i));
        }
    }
    
    model.setObjective(expr, GRB_MINIMIZE);
    model.optimize();
    cout<<"Min number of segments: "<< model.getObjective().getValue()<<endl;
    
    ofstream f_cut("min_segment_cuts.txt");
    for(int i=0;i<segments.size();i++){
        if(v_seg[i].get(GRB_DoubleAttr_X) == 1.0){
            f_cut << segments[i].p1.x<<' ' << segments[i].p1.y << ' ' << segments[i].p2.x <<' ' << segments[i].p2.y<<endl;
        }
    }
    f_cut.close();
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
    if(ptr!=ext_pts.begin() && !is_collision(segment(seg.p1, *(ptr-1)))) seg.p1 = *(--ptr);
    ptr = upper_bound(begin(ext_pts), end(ext_pts), seg.p2);
    if(ptr != ext_pts.end() && !is_collision(segment(seg.p2, *(ptr)))) seg.p2 = *ptr;
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

    polygon workspace;
    workspace.pts = vector<vec2>{vec2(l, d), vec2(r, d), vec2(r, u), vec2(l, u)};
    // obj.push_back(workspace);

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
                    // extend to both sides
                    extend(seg);
                    segments.push_back(seg);
                }
            }
            for(int k = i+1; k< obj.size(); k++){
                for(int l=0;l<obj[k].pts.size();l++){
                    segment seg(p, obj[k].pts[l]);
                    if(!is_collision(seg)){
                        extend(seg);
                        segments.push_back(seg);
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
        if(eq(segments[i].p1, segments[j].p1) && eq(segments[i].p2, segments[j].p2)) j++;
        else segments[++i] = segments[j++];
        if(j == segments.size()){
            segments.resize(i+1); break;
        }
    }
    for(const auto &seg: segments){
        f_seg << seg.p1.x <<' '<< seg.p1.y << ' '<<seg.p2.x <<' '<<seg.p2.y<<endl;
    }
    cell_construction();
    return 0;
}