# Marc Serrano Altena
# 29-11-2022
# this program simulates the gravitational pull of celestial body's

from matplotlib import pyplot as plt
from matplotlib import animation
from math import *
import random
from time import time

start_time = time()

class ball:
  dt = 0.01
  counter = 0

  particle_num = 0
  xmin_box = 0
  xmax_box = 10
  ymin_box = 0
  ymax_box = 10

  particles = []
  Nparticles = 3
  circles = []

  fig, ax = plt.subplots()
  ax.set_xlim([xmin_box, xmax_box])
  ax.set_ylim([ymin_box, ymax_box])
  ax.set_aspect('equal')


  def __init__(self):
    self.x = random.randint(ball.xmin_box + 1, ball.xmax_box - 1) 
    self.y = random.randint(ball.ymin_box + 1, ball.ymax_box - 1)
    
    self.r = random.uniform(0.05, 0.75)
    self.m = pi * (self.r * 3)**2 
    self.vx = random.uniform(-1.25, 1.25)
    self.vy = random.uniform(-1.25, 1.25)
    # the accelaration in the x and y direction
    self.a = random.uniform(0.5, 1)
    self.ax = 0
    self.ay = 0

    if self.particle_num == 0:
      self.plot = plt.Circle((self.x, self.y), radius = self.r, fc = 'red')

    else:
      self.plot = plt.Circle((self.x, self.y), radius = self.r, fc = 'blue')
    
    ball.ax.add_patch(self.plot)

    self.particle_num = ball.particle_num
    ball.particle_num += 1

  def __repr__(self):
    return f"ball_{self.particle_num}"
  

  # gives the speed decrease on collision
  def inelastic_col(self, direction):
    inelasticness = 0.25

    if direction == "x":

      if self.vx > inelasticness:
        self.vx -= inelasticness
      elif self.vx > 0:
        self.vx = 0

      elif self.vx < -inelasticness:
        self.vx += inelasticness
      elif self.vx < 0:
        self.vx = 0

    if direction == "y":

      if self.vy > inelasticness:
        self.vy -= inelasticness
      elif self.vy > 0:
        self.vy = 0

      elif self.vy < -inelasticness:
        self.vy += inelasticness
      elif self.vy < 0:
        self.vy = 0


  def update(self):

    # update atributes
    self.x += self.vx * ball.dt
    self.y += self.vy * ball.dt

    if not ((self.x <= ball.xmin_box + self.r + 0.01 and self.ax < 0) or (self.x > ball.xmax_box - self.r - 0.01 and self.ax > 0)):
      self.vx += self.ax * ball.dt

    if not ((self.y <= ball.xmin_box + self.r + 0.01 and self.ay < 0) or (self.y >= ball.xmax_box - self.r - 0.01 and self.ay > 0)):
      self.vy += self.ay * ball.dt

      
    # when hitting a wall; ± radius so that the bounce happens at the surface of the particle
    if self.x <= ball.xmin_box + self.r: 
        self.vx = - self.vx
        self.x += 0.01
        self.inelastic_col("x")
    elif self.x >= ball.xmax_box - self.r:
        self.vx = - self.vx
        self.x -= 0.01
        self.inelastic_col("x")

    if self.y <= ball.ymin_box + self.r: 
        self.vy = - self.vy
        self.y += 0.01
        self.inelastic_col("y")
    elif self.y >= ball.ymax_box - self.r:
        self.vy = - self.vy
        self.y -= 0.01
        self.inelastic_col("y")


    # for interactions with other particles
    for num2 in range(self.particle_num + 1, ball.Nparticles):
      # nx = (x2 - x1)
      nx = ball.particles[num2].x - self.x
      # ny = (y2 - y1)
      ny = ball.particles[num2].y - self.y
      # distance between center particles
      d = sqrt(nx**2 + ny**2)


      # gravitational pull between particles
      self.uax = nx / d
      self.uay = ny / d

      ball.particles[num2].uax = -self.uax
      ball.particles[num2].uay = -self.uay

      self.a = ball.particles[num2].m / d**2
      ball.particles[num2].a = self.m / d**2


      self.ax = self.a * self.uax
      self.ay = self.a * self.uay
      ball.particles[num2].ax = ball.particles[num2].a * ball.particles[num2].uax
      ball.particles[num2].ay = ball.particles[num2].a * ball.particles[num2].uay


      # collisions with other particles
      if d <= self.r + ball.particles[num2].r:

        # c is speed of particle 1, b the speed of particle 2
        cx = self.vx
        cy = self.vy
        bx = ball.particles[num2].vx
        by = ball.particles[num2].vy

        # m1 is mass of particle 1, m2 the mass of particle 2
        m1 = self.m
        m2 = ball.particles[num2].m

        # unit normal vector
        unx = nx/sqrt(nx**2 + ny**2)
        uny = ny/sqrt(nx**2 + ny**2)
        # unit tangent vector
        utx = -uny
        uty = unx

        # magnitude of the vectors before collision
        cn = cx*unx + cy*uny
        bn = bx*unx + by*uny

        ct = cx*utx + cy*uty
        bt = bx*utx + by*uty

        # magnitude vectors after collision (f --> final)
        cn_f = (cn*(m1 - m2) + 2 * m2 * bn) / (m1 + m2)
        bn_f = (bn*(m2 - m1) + 2 * m1 * cn) / (m1 + m2)
        ct_f = ct
        bt_f = bt

        # for particle 1
        self.vx = cn_f * unx + ct_f * utx
        self.vy = cn_f * uny + ct_f * uty

        # for particle 2
        ball.particles[num2].vx = bn_f * unx + bt_f * utx
        ball.particles[num2].vy = bn_f * uny + bt_f * uty


        # to avoid particles sticking together
        if nx > 0:
          ball.particles[num2].x += 0.005
          self.x -= 0.005

        else:
          self.x += 0.005
          ball.particles[num2].x -= 0.005

        if ny > 0:
          ball.particles[num2].y += 0.005
          self.y-= 0.005

        else:
          self.y += 0.005
          ball.particles[num2].y -= 0.005
        

          

  @staticmethod
  def animate(i_time):

    for num in range(ball.Nparticles):
      ball.particles[num].update()
      ball.circles[num].center = (ball.particles[num].x, ball.particles[num].y)


    return ball.circles



# main program
def main():

  for num in range(ball.Nparticles):
    ball.particles.append(ball())
    ball.circles.append(ball.particles[num].plot)
    

  #--/ do the actual animation (part of main program)
  anim = animation.FuncAnimation(ball.fig, ball.animate,
                              frames=20000, interval=2, blit=True)

  plt.show()



# main program
if __name__ == "__main__":
  main()