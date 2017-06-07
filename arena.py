import math
import ntpath
import random
import subprocess
import sys
import time

from plane import Plane, inted

class Player:
    def __init__(self, name):
        self.name = ntpath.basename(name)
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
        self.facing = (0,0)
        self.speed = [0,0]
        self.radius = 400
        self.friction_value = 0.85
        self.boost = 650
        self.has_boost = True
        self.shield_multiplier = 10
        self.shield_status = 0
        self.lost = False

    def update(self, commands):
        cmd = commands.split(" ")
        if cmd[2] == "BOOST":
            if self.has_boost == True: 
                cmd[2] = str(self.boost)
                self.has_boost == False
            else:
                self.lost = True
        elif cmd[2] == "SHIELD":
            if self.shield_status == 0:
                self.shield_status = 3
        if self.shield_status > 0:
            cmd[2] = 0
            self.shield_status -= 1
        c = [int(i) for i in cmd]
        world = Plane()
        world.moveTo((self.x,self.y))
        world.lookAt((self.facing[0],self.facing[1]))
        self.rotate(c, world)
        self.accelerate(c)
        self.move()
        self.friction()
        self.truncate()

    def rotate(self, cmd, world):
        world.lookAt18((cmd[0],cmd[1]))
        coords = world.realCoord((0,1))
        self.facing = float(coords[0][0]),float(coords[1][0])
        tan = -(float(coords[0][0])-self.x)/(float(coords[1][0])-self.y)
        self.angle = math.atan(tan)

    def accelerate(self, cmd):
        vecx, vecy = self.facing[0] - self.x, self.facing[1] - self.y
        vecx *= cmd[2]
        vecy *= cmd[2]
        self.speed[0] += vecx
        self.speed[1] += vecy

    def move(self):
        self.checkCollisions()
        self.x += self.speed[0]
        self.y += self.speed[1]

    def checkCollisions(self):
        pass

    def friction(self):
        self.speed[0] = self.speed[0] * self.friction_value
        self.speed[1] = self.speed[1] * self.friction_value

    def truncate(self):
        self.speed[0] = int(self.speed[0])
        self.speed[1] = int(self.speed[1])
        self.x = round(self.x)
        self.y = round(self.y)
        self.angle = int(self.angle)



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
        return (random.randrange(16000), random.randrange(9000))

    def initPods(self):
        self.Pods = [Pod() for i in range(4)]

    def start(self):
        self.sendPlayers("Laps: ", str(self.Laps))
        self.sendPlayers("Checkpoints: ", str(len(self.Checkpoints)))
        for cp in self.Checkpoints:
            self.sendPlayers("CP: ", str(cp[0]) + " " + str(cp[1]))
        world = Plane()
        world.moveTo(self.Checkpoints[0])
        world.lookAt(self.Checkpoints[1])
        for i in range(4):
            coord = inted(world.realCoord((-1500 + (1000*i),0)))
            self.Pods[i].x = coord[0][0]
            self.Pods[i].y = coord[1][0]
            self.Pods[i].facing = self.Checkpoints[1]
            atan_numerator = -(float(self.Checkpoints[1][0]) - coord[0][0])
            atan_denominator = (float(self.Checkpoints[1][1]) - coord[1][0])
            self.Pods[i].angle = math.atan(atan_numerator/atan_denominator)
        for i in range(10):
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
        msg1 = player.handle.stdout.readline().decode('utf-8').rstrip()
        print(player.name, msg1)
        msg2 = player.handle.stdout.readline().decode('utf-8').rstrip()
        print(player.name, msg2)
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
