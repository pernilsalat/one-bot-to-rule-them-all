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

