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

class Pod:
    def __init__(self):
        self.x = random.randrange(16000)
        self.y = random.randrange(9000)
        self.angle = 0.00
        self.facing = [0,0]
        self.speed = [0,0]
        self.radius = 400
        self.friction = 0.85
        self.boost = 650
        self.has_boost = True
        self.shield_multiplier = 10
        self.shield_status = 0
        self.lost = False

    def update(commands):
        cmd = commands.split(" ")
        if cmd[2] == "BOOST":
            if has_boost == True: 
                cmd[2] = str(self.boost)
                has_boost == False
            else:
                self.lost = True
        if cmd[2] == "SHIELD":
            if self.shield_status == 0:
                self.shield_status = 3
                cmd[2] = 0

class Game:
    def __init__(self, p1, p2):
        self.Player1 = Player(p1)
        self.Player2 = Player(p2)
        self.initCircuit()
        self.initPods()

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

    def initPods(self):
        self.Pods = [Pod() for i in range(4)]

    def start(self):
        self.sendPlayers("Laps: ", str(self.Laps))
        self.sendPlayers("Checkpoints: ", str(len(self.Checkpoints)))
        for cp in self.Checkpoints:
            self.sendPlayers("CP: ", str(cp[0]) + " " + str(cp[1]))
        self.sendTurn()
        self.calculate()

    def sendTurn(self):
        self.sendPod(0, self.Player1)
        self.sendPod(1, self.Player1)
        self.sendPod(2, self.Player1)
        self.sendPod(3, self.Player1)
        self.sendPod(2, self.Player2)
        self.sendPod(3, self.Player2)
        self.sendPod(0, self.Player2)
        self.sendPod(1, self.Player2)

    def calculate(self):
        bot0, bot1 = self.getResponse(self.Player1)
        bot2, bot3 = self.getResponse(self.Player2)
        self.Pods[0].update(bot0)
        self.Pods[0].update(bot1)
        self.Pods[0].update(bot2)
        self.Pods[0].update(bot3)

    def sendPlayers(self, debug, msg):
        self.sendPlayer1(debug,msg)
        self.sendPlayer2(debug,msg)
    
    def sendPlayer1(self, debug, msg):
        print(debug, msg)
        self.Player1.handle.stdin.write((msg+"\n").encode('utf-8'))

    def sendPlayer2(self, debug, msg):
        print(debug, msg)
        self.Player2.handle.stdin.write((msg+"\n").encode('utf-8'))

    def sendPod(self, index, player):
        pod = self.Pods[index]
        msg = str(pod.x) + " " + str(pod.y) + " "
        msg += str(pod.speed[0]) + " " + str(pod.speed[1]) + " "
        msg += str(int(pod.angle)) + " 1"
        print("Pod: ", msg)
        player.handle.stdin.write((msg+"\n").encode('utf-8'))

    def getResponse(self, player):
        player.handle.stdin.flush()
        msg1 = player.handle.stdout.readline()
        print(player.name, msg1.decode('utf-8'))
        msg2 = player.handle.stdout.readline()
        print(player.name, msg2.decode('utf-8'))
        return msg1, msg2


    def end(self):
        self.Player1.handle.kill()
        self.Player2.handle.kill()


if len(sys.argv) < 3:
    print("Usage: python3 arena.py Player1 Player2")
    sys.exit()


game = Game(sys.argv[1], sys.argv[2])
game.start()
game.end()
