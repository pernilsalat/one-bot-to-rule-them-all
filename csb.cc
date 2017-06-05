#pragma GCC optimize("Ofast", "inline", "omit-frame-pointer", "unroll-loops")

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

enum Type {CP, POD};

constexpr int DEPTH = 6;
constexpr float SHIELD_DEPTH = 9; //p is 1/S_D
constexpr float P_MUTATE_SHIELD = 15; //p for a single shield mutation. (*2 when considering both);
constexpr int POOL = 20;
constexpr int MAX_THRUST = 200;
constexpr int AMPLITUDE_ANGLE = 36;
constexpr int AMPLITUDE_MIN_THRUST = -30;
constexpr int AMPLITUDE_MAX_THRUST = 300;

constexpr float E = 0.00001;

float COS[360];
float SIN[360];

int iteration = -1;
int turn = 0;
int gen_ct = 0;
int best_gen_ct = 0;

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

inline int addrnd(int a, int b) {
    return a + rnd(b - a + 1);
}

//*****************************************************************************************//

class Collision {
    public:
        Unit* a;
        Unit* b;
        float time;

        Collision() {}

        Collision(Unit* a, Unit* b, float time) {
            this->a = a;
            this->b = b;
            this->time = time;
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
};

inline float dist(Point* p1, Point* p2) {
    return sqrt((p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y));
}

class Unit: public Point {
    private:
        float cache_x, cache_y, cache_vx, cache_vy;

    public:
        int id;
        Type type;
        float radius, vx, vy;

        inline virtual void bounce(Unit*){};

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
            float d = 1.0 / dist(this, p);

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
            this->radius = 600;
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

        int cache_ncpid, cache_checked, cache_shield;
        float cache_angle;
        bool cache_has_boost;

        Pod(int id) {
            this->id = id;
            this->radius = 400;
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
            
            vx += COS[(int)angle] * thrust;
            vy += SIN[(int)angle] * thrust;
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

            this->angle = angle;
            if (iteration == 0) this->angle = diff_angle(cps[1]);
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
                float score[ENEMIES];
                int thrusts0[2];
                int thrusts[2][ENEMIES][DEPTH-1];
                float angles0[2];
                float angles[2][ENEMIES][DEPTH-1];
                bool shieldTurn0[2];
                int shieldTurn[2][ENEMIES];

                Solution() {};

                inline void mutate () {
                        int r_turn = rnd(DEPTH);
                        randomize(r_turn);

                        int r = rnd(100);
                        if (r_turn) {
                                if (r < P_MUTATE_SHIELD) {
                                        int r_shield = addrnd(1, SHIELD_DEPTH);
                                        for (int i = 0; i < ENEMIES; ++i) shieldTurn[0][i] = r_shield;
                                } else if (r < P_MUTATE_SHIELD * 2) {
                                        int r_shield = addrnd(1, SHIELD_DEPTH);
                                        for (int i = 0; i < ENEMIES; ++i) shieldTurn[1][i] = r_shield;
                                }
                        } else {
                                shieldTurn0[0] = rnd(SHIELD_DEPTH) == 0;
                                shieldTurn0[1] = rnd(SHIELD_DEPTH) == 0;
                        }

                        for (int i = 0; i < ENEMIES; ++i) score[i] = -1;
                }

                inline void randomize (int idx, bool full = false) {
                        if (idx) {
                        } else {
                                if (full || r == 0) {
                                        angles0[0] = max(-18, min(18, addrnd(-AMPLITUDE_ANGLE, AMPLITUDE_ANGLE)));
                                        angles0[1] = max(-18, min(18, addrnd(-AMPLITUDE_ANGLE, AMPLITUDE_ANGLE)));
                                }
                                if (full || r == 1) {
                                        thrusts0[0] = max(0, min(MAX_THRUST, addrnd(AMPLITUDE_MIN_THRUST, AMPLITUDE_MAX_THRUST)));
                                        thrusts0[1] = max(0, min(MAX_THRUST, addrnd(AMPLITUDE_MIN_THRUST, AMPLITUDE_MAX_THRUST)));
                                }
                        }
                }

                inline void randomize () {
                        for (int i = 0; i < DEPTH; ++i) randomize(i, TRUE);
                        for (int i = 0; i < ENEMIES; ++i) {
                                shieldTurn[0][i] = addrnd(1, SHIELD_DEPTH);
                                shieldTurn[1][i] = addrnd(1, SHIELD_DEPTH);
                        }
                        shieldTurn0[0] = rnd(SHIELD_DEPTH) == 0;
                        shieldTurn0[1] = rnd(SHIELD_DEPTH) == 0;

                        for (int i = 0; i < ENEMIES; ++i) score[i] = -1;
                }
};

//*****************************************************************************************//

class Bot {
    public:
        int id = 0;

        Bot() {};

        Bot(int id) {
            this->id = id;
        }

        inline virtual void move(){};

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

            int thrust = fabs(raw_angle) < 60 ? MAX_THRUST : 100;
            thrust = fabs(raw_angle) < 90 ? thrust : 0;
            float angle = max((float) -18, min((float) 18, raw_angle));

            if (!for_output) pod->apply(thrust, angle);
            else print_move(-1, thrust, angle, pod);
        }

        inline void move_blocker(bool for_output = false) {
            Pod* pod = !for_output ? blocker() : pods[1];
            // hisRunner is not used, is this a mistake?
            // Pod* hisRunner = runner(2-id);

            Checkpoint* cp = cps[pod->ncpid];
            Point t(cp->x - 3*pod->vx, cp->y - 3*pod->vy);
            float raw_angle = pod->diff_angle(&t);

            int thrust = fabs(raw_angle) < 60 ? MAX_THRUST : 100;
            thrust = fabs(raw_angle) < 90 ? thrust : 0;
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
                if (col_time > -1 && col_time + t < 1.0 && (first_col.time == -1 || col_time < first_col.time)) {
                    first_col.a = pods[i];
                    first_col.b = pods[j];
                    first_col.time = col_time;
                }
            }

        }

        if (first_col.time == -1) {
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
                if (col_time > t && col_time < first_col.time) {
                    pods[i]->bounce(cps[pods[i]->ncpid]);
                }

                pods[i]->move(first_col.time);
            }

            first_col.a->bounce(first_col.b);
            t += first_col.time;
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
    
    for (int i = 0; i < 360; ++i) {
        float ra = i * M_PI / 180.0;
        COS[i] = cos(ra);
        SIN[i] = sin(ra);
    }
    
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

    SearchBot me;
    me.oppBots.push_back(&opp_reflex);

    while (1) {
        iteration++;


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

        float time_limit = iteration ? 0.142 : 0.98;

        save();

        //TODO Generate opponents

        me.solve(time_limit, iteration > 0);

        if (iteration) {
            cerr << "Avg. simulations: " << gen_ct * POOL * DEPTH / iteration << endl;
            cerr << "Avg. generations: " << gen_ct / iteration << endl;
            cerr << "Best generation:  " << best_gen_ct / iteration << endl;
        } else {
            gen_ct = 0;
            best_gen_ct = 0;
        }

        print_move(me.sol.shieldTurn1, me.sol.thrusts[0], me.sol.angles[0], pods[0]);
        print_move(me.sol.shieldTurn2, me.sol.thrusts[DEPTH], me.sol.angles[DEPTH], pods[1]);

    }
}
