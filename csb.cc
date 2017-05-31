#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("omit-frame-pointer")
#pragma GCC optimize("unroll-loops")

#include "stdio.h"
#include "math.h"
#include <iostream>
#include <algorithm>
#include <memory>
#include <chrono>
#include <vector>

using namespace std;
using namespace std::chrono;

high_resolution_clock::time_point now = high_resolution_clock::now();
#define TIME duration_cast<duration<double>>(high_resolution_clock::now() - now).count()

class Point;
class Unit;
class Pod;
class Collision;
class Checkpoint;
class Solution;
class Bot;
void save();
void load();
void play();
void print_move(int, int, float, Pod*);

constexpr int CP  = 0;
constexpr int POD = 1;

constexpr int DEPTH = 6;
constexpr float SHIELD_DEPTH = 9; //p is 1/S_D
constexpr float P_MUTATE_SHIELD = 15; //p for a single shield mutation. (*2 when considering both);
constexpr int POOL = 20;
constexpr int MAX_THRUST = 200;

constexpr float E = 0.00001;

int r = -1;
int turn = 0;
int gen_ct = 0;
int best_gen_ct = 0;
bool is_p2 = false;

int cp_ct, laps;
int myTimeout = 100, hisTimeout = 100;
int cache_myTimeout, cache_hisTimeout;

Pod* pods[4];
Checkpoint* cps[10];

//*****************************************************************************************//

inline int fastrand() {
    static unsigned int g_seed = 42;
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

inline int rnd(int b) {
    return fastrand() % b;
}

inline int rnd(int a, int b) {
    return a + rnd(b - a + 1);
}

//*****************************************************************************************//

class Collision {
public:
    Unit* a;
    Unit* b;
    float t;

    Collision() {}

    Collision(Unit* a, Unit* b, float t) {
        this->a = a;
        this->b = b;
        this->t = t;
    }
};

//*****************************************************************************************//

class Point {
public:
    float x, y;

    Point() {};

    Point(float x, float y) {
        this->x = x;
        this->y = y;
    }

    inline Point closest(Point* a, Point* b) {
        float da = b->y - a->y;
        float db = a->x - b->x;
        float c1 = da*a->x + db*a->y;
        float c2 = -db*x + da*y;
        float det = da*da + db*db;

        float cx, cy;
        if (det != 0) {
            cx = (da*c1 - db*c2) / det;
            cy = (da*c2 + db*c1) / det;
        } else {
            cx = x, cy = y;
        }

        return Point(cx, cy);
    }
};

inline float dist(float x1, float y1, float x2, float y2) {
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

inline float dist2(float x1, float y1, float x2, float y2) {
	return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
}

inline float dist(Point* p, float x, float y) {
	return sqrt((x - p->x) * (x - p->x) + (y - p->y) * (y - p->y));
}

inline float dist2(Point* p, float x, float y) {
	return (x - p->x) * (x - p->x) + (y - p->y) * (y - p->y);
}

inline float dist(Point* p1, Point* p2) {
	return sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
}

inline float dist2(Point* p1, Point* p2) {
	return (p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y);
}

class Unit: public Point {
private:
    float cache_x, cache_y, cache_vx, cache_vy;

public:
    int id, type;
    float r, vx, vy;

    inline virtual void bounce(Unit* u) {};

    inline float collision_time(Unit* u) {
        if (vx == u->vx && vy == u->vy) {
            return -1;
        }

        float sr2 = u->type == CP ? 357604 : 640000;

        float dx = x - u->x;
        float dy = y - u->y;
        float dvx = vx - u->vx;
        float dvy = vy - u->vy;
        float a = dvx*dvx + dvy*dvy;

        if (a < E) return -1;

        float b = -2.0*(dx*dvx + dy*dvy);
        float delta = b*b - 4.0*a*(dx*dx + dy*dy - sr2);

        if (delta < 0.0) return -1;

        float t = (b - sqrt(delta))*(1.0/(2.0*a));

        if (t <= 0.0 || t > 1.0) return -1;

        return t;
    }
    
    float speedTo(Point* p) {
        float d = 1.0 / dist(this, p->x, p->y);

        float dx = (p->x - this->x) * d;
        float dy = (p->y - this->y) * d;
        float nspeed = vx*dx + vy*dy;
        float ospeed = dy*vx - dx*vy;

        return nspeed - (5 * ospeed * ospeed * d); //Magus used nspeed - (5 * ospeed * ospeed * d)
    }
    
    inline void save() {
        cache_x = x;
        cache_y = y;
        cache_vx = vx;
        cache_vy = vy;
    }

    inline void load() {
        x = cache_x;
        y = cache_y;
        vx = cache_vx;
        vy = cache_vy;
    }
};

//*****************************************************************************************//

class Checkpoint: public Unit {
public:
    Checkpoint(int id, float x, float y) {
        this->id = id;
        this->x = x;
        this->y = y;

        this->vx = this->vy = 0;
        this->type = CP;
        this->r = 600;
    }

    inline void bounce(Unit*) {}
};

class Pod: public Unit {
public:
    float angle = -1;
    float next_angle = -1;
    bool has_boost;
    int ncpid, checked, shield;
    Pod* partner;

    float cache_ncpid, cache_checked, cache_timeout, cache_shield, cache_angle, cache_has_boost;

    Pod(int id) {
        this->id = id;
        this->r = 400;
        this->type = POD;
        this->ncpid = 1;
        this->has_boost = true;
        this->checked = this->shield = 0;
    }
    
    inline float score() {    
        return checked*500000 - dist(this, cps[this->ncpid]) + this->speedTo(cps[this->ncpid]);
    }

    inline void apply(int thrust, float angle) {
		
        this->angle += angle;
        if (this->angle >= 360.) {
            this->angle = this->angle - 360.;
        } else if (this->angle < 0.0) {
            this->angle += 360.;
        }

        if (thrust == -1) {
            this->shield = 4;
        } else {
            boost(thrust);
        }
    }

    inline void rotate(Point* p) {
        float a = diff_angle(p);
        a = max((float)-18., min((float)18., a));

        angle += a;
        if (angle >= 360.) {
            angle = angle - 360.;
        } else if (angle < 0.0) {
            angle += 360.;
        }
    }

    inline void boost(int thrust) {
        if (shield > 0) return;

        float ra = angle * M_PI / 180.0;

        vx += cos(ra) * thrust;
        vy += sin(ra) * thrust;
    }

    inline void move(float t) {
        x += vx * t;
        y += vy * t;
    }

    inline void end() {
        x = round(x);
        y = round(y);
        vx = trunc(vx * 0.85);
        vy = trunc(vy * 0.85);

        if (checked >= cp_ct * laps) {
            ncpid = 0;
            checked = cp_ct * laps;
        }

        if (shield > 0) shield--;
    }

    inline void bounce(Unit* u) {
        if (u->type == CP) {
            checked += 1;
            if (id == 0 || id == 1) { //I put at 101 because I will then take one out.
				myTimeout = 101;
			} else {
				hisTimeout = 101;
			}
            ncpid = (ncpid + 1 == cp_ct) ? 0 : ncpid + 1;
            return;
        }

        bounce_w_pod(static_cast<Pod*>(u));
    }

    inline void bounce_w_pod(Pod* u) {
        float m1 = shield == 4 ? 10. : 1.;
        float m2 = u->shield == 4 ? 10. : 1.;
        float mcoeff = (m1 + m2) / (m1 * m2);

        float nx = x - u->x;
        float ny = y - u->y;
        float dst2 = nx*nx + ny*ny;
        float dvx = vx - u->vx;
        float dvy = vy - u->vy;
        float prod = (nx*dvx + ny*dvy) / (dst2 * mcoeff);
        float fx = nx * prod;
        float fy = ny * prod;
        float m1_inv = 1.0 / m1;
        float m2_inv = 1.0 / m2;

        vx -= fx * m1_inv;
        vy -= fy * m1_inv;
        u->vx += fx * m2_inv;
        u->vy += fy * m2_inv;

        float impulse = sqrt(fx*fx + fy*fy);
        if (impulse < 120.) {
            float df = 120.0 / impulse;
            fx *= df;
            fy *= df;
        }

        vx -= fx * m1_inv;
        vy -= fy * m1_inv;
        u->vx += fx * m2_inv;
        u->vy += fy * m2_inv;
    }

    inline float diff_angle(Point* p) {
        float a = get_angle(p);
        float right = angle <= a ? a - angle : 360. - angle + a;
        float left = angle >= a ? angle - a : angle + 360. - a;

        if (right < left) {
            return right;
        }

        return -left;
    }

    inline float get_angle(Point* p) {
        float d = dist(this, p);
        float dx = (p->x - x) / d;
        float dy = (p->y - y) / d;

        float a = acos(dx) * 180 / M_PI;

        if (dy < 0) {
            a = 360 - a;
        }

        return a;
    }

    inline void update(int x, int y, int vx, int vy, float angle, int ncpid) {
        if (shield > 0) shield--;
        if (ncpid != this->ncpid) {
            ++checked;
        } 

        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
        this->ncpid = ncpid;

        if (is_p2 && id > 1) swap(angle, this->next_angle);
        this->angle = angle;
        if (::r == 0) this->angle = diff_angle(cps[1]);
    }

    inline void update(int shield, bool has_boost) {
        this->shield = shield;
        this->has_boost = has_boost;
    }

    inline void save() {
        Unit::save();
        cache_ncpid = ncpid;
        cache_checked = checked;
        cache_shield = shield;
        cache_angle = angle;
        cache_has_boost = has_boost;
    }

    inline void load() {
        Unit::load();
        ncpid   = cache_ncpid;
        checked = cache_checked;
        shield  = cache_shield;
        angle   = cache_angle;
        has_boost = cache_has_boost;
    }
};

//*****************************************************************************************//

class Solution {
public:
    float score = -1;
    int thrusts[DEPTH*2];
    float angles[DEPTH*2];
	int shieldTurn1, shieldTurn2;

    Solution(bool with_rnd = false) {
        if (with_rnd) randomize();
    }
    
    inline void shift() {
        for (int i = 1; i < DEPTH; ++i) {
            angles[i-1]        = angles[i];
            thrusts[i-1]       = thrusts[i];
            angles[i-1+DEPTH]  = angles[i+DEPTH];
            thrusts[i-1+DEPTH] = thrusts[i+DEPTH];
        }
        randomize(DEPTH-1, true);
		
		shieldTurn1 = (shieldTurn1 == 0) ? rnd(SHIELD_DEPTH) : shieldTurn1 - 1;
		shieldTurn2 = (shieldTurn2 == 0) ? rnd(SHIELD_DEPTH) : shieldTurn2 - 1;
		
        score = -1;
    }

    inline void mutate() {
        randomize(rnd(DEPTH));
		int r = rnd(100);
		if (r < P_MUTATE_SHIELD) {
			shieldTurn1 = rnd(SHIELD_DEPTH);
		} else if (r < P_MUTATE_SHIELD * 2) {
			shieldTurn2 = rnd(SHIELD_DEPTH);
		}
    }

    inline void randomize(int idx, bool full = false) {
        int r = rnd(2);
        if (full || r == 0) {
            angles[idx] = max(-18, min(18, rnd(-22, 22)));
            angles[idx+DEPTH] = max(-18, min(18, rnd(-22, 22)));
        }

        if (full || r == 1) {
            thrusts[idx] = max(0, min(MAX_THRUST, rnd((int) (-0.1*MAX_THRUST), 1.5*MAX_THRUST)));
            thrusts[idx+DEPTH] = max(0, min(MAX_THRUST, rnd((int) (-0.1*MAX_THRUST), 1.5*MAX_THRUST)));
        }
        score = -1;
    }

    inline void randomize() {
        for (int i = 0; i < DEPTH; ++i) randomize(i, true);
		shieldTurn1 = rnd(SHIELD_DEPTH);
		shieldTurn2 = rnd(SHIELD_DEPTH);
    }
};

inline Solution merge(Solution* a, Solution* b) {
    Solution child = a;
    
    for (int i = 0; i < DEPTH; i++) {
        if (rnd(2) == 0) {
            child.angles[i] = b->angles[i];
            child.angles[i+DEPTH] = b->angles[i+DEPTH];
            child.thrusts[i] = b->thrusts[i];
            child.thrusts[i+DEPTH] = b->thrusts[i+DEPTH];
        }
    }
    
    switch (rnd(4)) {
        case 0:
            child.shieldTurn1 = b->shieldTurn1;
            break;
        case 1:
            child.shieldTurn2 = b->shieldTurn2;
            break;
        case 2:
            child.shieldTurn1 = b->shieldTurn1;
            child.shieldTurn2 = b->shieldTurn2;
            break;
    }
    
    child.score = -1;
    
    return child;
}

//*****************************************************************************************//

class Bot {
public:
    int id = 0;

    Bot() {};

    Bot(int id) {
        this->id = id;
    }

    inline virtual void move() = 0;

    inline Pod* runner() {
        return runner(pods[id], pods[id+1]);
    }

    inline Pod* blocker() {
        return blocker(pods[id], pods[id+1]);
    }
    
    inline Pod* runner(int i) {
        return runner(pods[i], pods[i+1]);
    }

    inline Pod* blocker(int i) {
        return blocker(pods[i], pods[i+1]);
    }

    inline Pod* runner(Pod* pod0, Pod* pod1) {
        if (pod0->checked > pod1->checked) return pod0;
        if (pod1->checked > pod0->checked) return pod1;
        return pod0->score() - pod1->score() >= -1000 ? pod0 : pod1;
    }

    inline Pod* blocker(Pod* pod0, Pod* pod1) {
        return runner(pod0, pod1)->partner;
    }
};

class ReflexBot : public Bot {
public:
    ReflexBot() {}

    ReflexBot(int id) {
        this->id = id;
    }

    inline void move() {
        move_runner();
        move_blocker();
    }

    inline void move_as_main() {
        move_runner(true);
        move_blocker(true);
    }

    inline void move_runner(bool for_output = false) {
        Pod* pod = !for_output ? runner() : pods[0];

        Checkpoint* cp = cps[pod->ncpid];
        Point t(cp->x - 3*pod->vx, cp->y - 3*pod->vy);
        float raw_angle = pod->diff_angle(&t);

        int thrust = abs(raw_angle) < 60 ? MAX_THRUST : 100;
        thrust = abs(raw_angle) < 90 ? thrust : 0;
        float angle = max((float) -18, min((float) 18, raw_angle));

        if (!for_output) pod->apply(thrust, angle);
        else print_move(-1, thrust, angle, pod);
    }

    inline void move_blocker(bool for_output = false) {
        Pod* pod = !for_output ? blocker() : pods[1];
        Pod* hisRunner = runner(2-id);

        Checkpoint* cp = cps[pod->ncpid];
        Point t(cp->x - 3*pod->vx, cp->y - 3*pod->vy);
        float raw_angle = pod->diff_angle(&t);

        int thrust = abs(raw_angle) < 60 ? MAX_THRUST : 100;
        thrust = abs(raw_angle) < 90 ? thrust : 0;
        float angle = max((float) -18, min((float) 18, raw_angle));

        if (!for_output) pod->apply(thrust, angle);
        else print_move(-1, thrust, angle, pod);
    }
};

class SearchBot : public Bot {
public:
    Solution sol;
    vector<Bot*> oppBots;

    SearchBot() {}

    SearchBot(int id) {
        this->id = id;
    }

    inline void move(Solution* sol) {
		if (sol->shieldTurn1 == 0) {
			pods[id]->apply(-1, sol->angles[turn]);
		} else {
			pods[id]->apply(sol->thrusts[turn], sol->angles[turn]);
		}
		if (sol->shieldTurn2 == 0) {
			pods[id+1]->apply(-1, sol->angles[turn+DEPTH]);
		} else {
			pods[id+1]->apply(sol->thrusts[turn+DEPTH], sol->angles[turn+DEPTH]);
		}
    }

    inline void move() {
        move(&sol);
    }

    inline void solve(float time, bool with_seed = false) {
        Solution best;
        Solution* pool = new Solution[POOL];
        Solution* newPool = new Solution[POOL];
        int pool_ct = 0;
        
        int current_gen = 0;
        int best_gen = 1;
        
        if (with_seed) {
            best = sol;
            best.shift();
            pool[pool_ct++] = best;
            while (pool_ct < POOL/5) {
                Solution bestProxy = sol;
                bestProxy.shift();
                pool[pool_ct++] = bestProxy;
                
                if (get_score(&bestProxy) > get_score(&best)) best = bestProxy;
            }
        } else {
            best.randomize();
            if (r == 0 && dist(pods[id], cps[1]) > 4000) best.thrusts[0] = 650;
            pool[pool_ct++] = best;
        }
            
        while (pool_ct < POOL) {
            Solution rSol;
            rSol.randomize();
            if (get_score(&rSol) > get_score(&best)) best = rSol;
            pool[pool_ct++] = rSol;
        }

        ++gen_ct;
        ++current_gen;

        //start creating new generations
        while (TIME < time) {
            newPool[0] = best;
            Solution mBest = best;
            mBest.mutate();
            if (get_score(&mBest) > get_score(&best)) best = mBest;
            newPool[1] = mBest;
            
            int poolC = 2;
            
            while (TIME < time && poolC < POOL) {
                int aIdx = rnd(POOL);
                int bIdx = rnd(POOL);
                while (aIdx == bIdx) bIdx = rnd(POOL);
                
                int firstIdx = (get_score(&pool[aIdx]) > get_score(&pool[bIdx])) ? aIdx : bIdx;
                
                aIdx = rnd(POOL);
                while (aIdx == firstIdx) aIdx = rnd(POOL); 
                bIdx = rnd(POOL);
                while (aIdx == bIdx || bIdx == firstIdx) bIdx = rnd(POOL);
                
                int secondIdx = (get_score(&pool[aIdx]) > get_score(&pool[bIdx])) ? aIdx : bIdx;
                
                Solution child = merge(&pool[firstIdx], &pool[secondIdx]);
                child.mutate();
                
                if (get_score(&child) > get_score(&best)) {
                    best = child;
                    best_gen = current_gen;
                }
                
                newPool[poolC++] = child;
            }
            
            for (int i = 0; i < POOL; ++i) pool[i] = newPool[i];
            
            ++gen_ct;
            ++current_gen;
        }
        
        best_gen_ct += best_gen;
        sol = best;
    }

    inline float get_score(Solution* sol) {
        if (sol->score == -1) {
            vector<float> scores;
            for (Bot* oppBot : oppBots) {
                scores.push_back(get_bot_score(sol, oppBot));
            }

            sol->score = *min_element(scores.begin(), scores.end());
        }

        return sol->score;
    }

    inline float get_bot_score(Solution* sol, Bot* opp) {
        float score = 0;
        while (turn < DEPTH) {
            move(sol);
            opp->move();
            play();
            if (turn == 0) score += 0.1*evaluate();
            ++turn;
        }
        score += 0.9*evaluate();
        load();

        return score;
    }

    inline float evaluate() {
        Pod* my_runner = runner(pods[id], pods[id+1]);
        Pod* my_blocker = blocker(pods[id], pods[id+1]);
        Pod* opp_runner = (id == 0) ? runner(pods[(id+2)], pods[(id+3)]) : runner(pods[0], pods[1]);
        Pod* opp_blocker = (id == 0) ? blocker(pods[(id+2)], pods[(id+3)]) : blocker(pods[0], pods[1]);

        if (myTimeout <= 0) return -1e7;
        if (hisTimeout <= 0) return 1e7;
        if (opp_runner->checked == laps*cp_ct || opp_blocker->checked == laps*cp_ct) return -1e7;
        if (my_runner->checked == laps*cp_ct || my_blocker->checked == laps*cp_ct) return 1e7;

        float score = my_runner->score() - opp_runner->score();
                
        float d_or = dist(opp_runner, cps[opp_runner->ncpid]);
        float d_mb = dist(my_blocker, cps[opp_runner->ncpid]);
        
        if (d_or < d_mb) {
            int nncpid = (opp_runner->ncpid == cp_ct - 1) ? 0 : opp_runner->ncpid + 1;
            score -= dist(my_blocker, cps[nncpid]) * 0.5;
            score -= 500;
        } else {
            score -= dist(my_blocker, cps[opp_runner->ncpid]) * 0.5;
        }
        score -= abs(my_blocker->diff_angle(opp_runner)) * 0.5;
        
        //check angl between opp blocker and my next cp
        float x1 = cps[my_runner->ncpid]->x - my_runner->x;
        float y1 = cps[my_runner->ncpid]->y - my_runner->y;
        float x2 = opp_blocker->x - my_runner->x;
        float y2 = opp_blocker->y - my_runner->y;
        
        float dot = x1*x2 + y1*y2;      // dot product
        float det = x1*y2 - y1*x2;      // determinant
        float a_bcp = atan2(det, dot); // atan2(y, x) or atan2(sin, cos)
        
        if (a_bcp > M_PI) a_bcp = 2 * M_PI - a_bcp;
        
        score += a_bcp * 5;
        score -= my_runner->shield * 10;
        
        //now for opp
        x1 = cps[opp_runner->ncpid]->x - opp_runner->x;
        y1 = cps[opp_runner->ncpid]->y - opp_runner->y;
        x2 = my_blocker->x - opp_runner->x;
        y2 = my_blocker->y - opp_runner->y;
        
        dot = x1*x2 + y1*y2;      // dot product
        det = x1*y2 - y1*x2;      // determinant
        a_bcp = atan2(det, dot); // atan2(y, x) or atan2(sin, cos)
        
        if (a_bcp > M_PI) a_bcp = 2 * M_PI - a_bcp;
        
        score -= a_bcp * 5;
        
        
        return score;
    }
};

//*****************************************************************************************//

inline void save() {
	for (int i = 0; i < 4; ++i) pods[i]->save();
	cache_myTimeout = myTimeout;
	cache_hisTimeout = hisTimeout;
}

inline void load() {
    for (int i = 0; i < 4; ++i) pods[i]->load();
    turn = 0;
	myTimeout = cache_myTimeout;
	hisTimeout = cache_hisTimeout;
}

inline void play() {
    float t = 0.0;
    while (t < 1.0) {
        Collision first_col = {NULL, NULL, -1};
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                float col_time = pods[i]->collision_time(pods[j]);
                if (col_time > -1 && col_time + t < 1.0 && (first_col.t == -1 || col_time < first_col.t)) {
                    first_col.a = pods[i];
                    first_col.b = pods[j];
                    first_col.t = col_time;
                }
            }

        }

        if (first_col.t == -1) {
            for (int i = 0; i < 4; ++i) {
                
                float col_time = pods[i]->collision_time(cps[pods[i]->ncpid]);
                if (col_time > t && col_time + t < 1.0) {
                    pods[i]->bounce(cps[pods[i]->ncpid]);
                }
                
                
                pods[i]->move(1.0 - t);
            }
            t = 1.0;
        } else {
            for (int i = 0; i < 4; ++i) {
               
               float col_time = pods[i]->collision_time(cps[pods[i]->ncpid]);
                if (col_time > t && col_time < first_col.t) {
                    pods[i]->bounce(cps[pods[i]->ncpid]);
                }
                
                pods[i]->move(first_col.t);
            }

            first_col.a->bounce(first_col.b);
            t += first_col.t;
        }
    }

	--myTimeout;
	--hisTimeout;
	
    for (int i = 0; i < 4; ++i) {
        pods[i]->end();
    }
}

inline void print_move(int shieldTurn, int thrust, float angle, Pod* pod) {
    float a = pod->angle + angle;

    if (a >= 360.0) {
        a = a - 360.0;
    } else if (a < 0.0) {
        a += 360.0;
    }

    a = a * M_PI / 180.0;
    float px = pod->x + cos(a) * 10000.0;
    float py = pod->y + sin(a) * 10000.0;

    if (shieldTurn == 0) {
        printf("%d %d SHIELD\n", (int) round(px), (int) round(py));
        pod->shield = 4;
    } else if (thrust == 650) {
        pod->has_boost = false;
        printf("%d %d BOOST\n", (int) round(px), (int) round(py));
    } else {
        printf("%d %d %d\n", (int) round(px), (int) round(py), thrust);
    }
}

//*****************************************************************************************//

int main() {
    string str1 = "";
    string str2 = "";
    
    cin >> laps >> cp_ct;
    for (int i = 0; i < cp_ct; i++) {
        int cx, cy;
        cin >> cx >> cy;
        cps[i] = new Checkpoint(i, cx, cy);
    }

    for (int i = 0; i < 4; i++) pods[i] = new Pod(i);

    pods[0]->partner = pods[1];
    pods[1]->partner = pods[0];
    pods[2]->partner = pods[3];
    pods[3]->partner = pods[2];

    ReflexBot me_reflex;
    ReflexBot opp_reflex(2);

    SearchBot opp(2);
    opp.oppBots.push_back(&me_reflex);

    SearchBot me;
    me.oppBots.push_back(&opp);
    me.oppBots.push_back(&opp_reflex);

    while (1) {
        r++;

        
        int x1, y1, vx1, vy1, angle1, ncpid1;
        int x2, y2, vx2, vy2, angle2, ncpid2;
        int x3, y3, vx3, vy3, angle3, ncpid3;
        int x4, y4, vx4, vy4, angle4, ncpid4;
        cin >> x1 >> y1 >> vx1 >> vy1 >> angle1 >> ncpid1;
        cin >> x2 >> y2 >> vx2 >> vy2 >> angle2 >> ncpid2;
        cin >> x3 >> y3 >> vx3 >> vy3 >> angle3 >> ncpid3;
        cin >> x4 >> y4 >> vx4 >> vy4 >> angle4 >> ncpid4;
        if (ncpid3 == -1 || ncpid4 == -1) while(1) printf("8000 4500 0\n");
            
        if (pods[0]->ncpid != ncpid1 || pods[1]->ncpid != ncpid2) {
            myTimeout = 100;
        } else {
            --myTimeout;
        }
        
        if (pods[2]->ncpid != ncpid3 || pods[3]->ncpid != ncpid4) {
            hisTimeout = 100;
        } else {
            --hisTimeout;
        }     
            
        pods[0]->update(x1, y1, vx1, vy1, angle1, ncpid1);
        pods[1]->update(x2, y2, vx2, vy2, angle2, ncpid2);
        pods[2]->update(x3, y3, vx3, vy3, angle3, ncpid3);
        pods[3]->update(x4, y4, vx4, vy4, angle4, ncpid4);
        

        now = high_resolution_clock::now();

        float time_limit = r ? 0.142 : 0.98;
		
		save();

        // use this to test reflex bot behavior
        //me_reflex.move_as_main();

        opp.solve(time_limit*0.2, r > 0);
        me.solve(time_limit, r > 0);
        
        if (r) {
            cerr << "Avg. simulations: " << gen_ct * POOL * DEPTH / r << endl;
            cerr << "Avg. generations: " << gen_ct / r << endl;
            cerr << "Best generation:  " << best_gen_ct / r << endl;
        } else {
            gen_ct = 0;
            best_gen_ct = 0;
        }
    
        str1 = str1 + "(" + to_string(me.sol.thrusts[0]) + "," + to_string((int)me.sol.angles[0]) + ")";
        str2 = str2 + "(" + to_string(me.sol.thrusts[1]) + "," + to_string((int)me.sol.angles[1]) + ")";
        
        cerr << str1 << endl;
        cerr << str2 << endl;
        
        print_move(me.sol.shieldTurn1, me.sol.thrusts[0], me.sol.angles[0], pods[0]);
        print_move(me.sol.shieldTurn2, me.sol.thrusts[DEPTH], me.sol.angles[DEPTH], pods[1]);

    }
}