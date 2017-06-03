import subprocess
import sys
import time
import random

class Player:
    def __init__(self, name):
        self.name = name
        self.handle = subprocess.Popen(name,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    def setpos(self, x, y):
        self.x = x
        self.y = y

class Game:
    def __init__(self, p1, p2):
        self.Player1 = Player(p1)
        self.Player2 = Player(p2)
        self.initCircuit()

    def initCircuit(self):
        self.initRand()
        self.initCheckpoints()

    def initRand(self):
        self.seed = time.time()
        print("Seed is:",self.seed)
        random.seed(self.seed)

    def initCheckpoints(self):
        N = random.randrange(2, 9)
        self.Laps = random.randrange(2,7)
        self.Checkpoints = [self.initCheckpoint() for i in range(N)]

    def initCheckpoint(self):
        return (random.randrange(16000), random.randrange(9000), 600)

    def start(self):
        print("Laps: ", self.Laps)
        print("Checkpoints: ", len(self.Checkpoints))
        for cp in self.Checkpoints:
            print(cp[0],cp[1])
        pass

    def end(self):
        self.Player1.handle.kill()
        self.Player2.handle.kill()


if len(sys.argv) < 3:
    print("Usage: python3 arena.py Player1 Player2")
    sys.exit()


game = Game(sys.argv[1], sys.argv[2])
game.start()
game.end()
