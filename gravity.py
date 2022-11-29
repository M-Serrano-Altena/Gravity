# Marc Serrano Altena
# 29-11-2022
# this program simulates the gravitational pull of celestial body's

from matplotlib import pyplot as plt
from matplotlib import animation
from math import *
import random

class ball:
  counter = 0
  particle_num = 0
  xmin_box = 0
  xmax_box = 10
  ymin_box = 0
  ymax_box = 10

  particles = []
  Nparticles = 20
  circles = []

  fig, ax = plt.subplots()
  ax.set_xlim([xmin_box, xmax_box])
  ax.set_ylim([ymin_box, ymax_box])
  ax.set_aspect('equal')


  def __init__(self):
    self.x = random.randint(ball.xmin_box + 1, ball.xmax_box - 1) 
    self.y = random.randint(ball.ymin_box + 1, ball.ymax_box - 1)
    # random angle the balls start at
    self.a = random.uniform(0, 2 * pi)
    self.v = random.uniform(0.5, 1.5)
    self.r = random.uniform(0.15, 0.5)
    self.m = pi * self.r**2
    self.vx = self.v*cos(self.a)
    self.vy = self.v*sin(self.a)

    if self.particle_num == 0:
      self.plot = plt.Circle((self.x, self.y), radius = self.r, fc = 'red')

    else:
      self.plot = plt.Circle((self.x, self.y), radius = self.r, fc = 'blue')
    
    ball.ax.add_patch(self.plot)

    self.particle_num = ball.particle_num
    ball.particle_num += 1

  def __repr__(self):
    return f"ball_{self.particle_num}"


  def update(self):
    dt = 0.01

    # update atributes
    self.x += self.vx * dt
    self.y += self.vy * dt

    # when hitting a wall; ± radius so that the bounce happens at the surface of the particle
    if self.x <= ball.xmin_box + self.r: 
        self.vx = - self.vx
        self.x += 0.01
    elif self.x >= ball.xmax_box - self.r:
        self.vx = - self.vx
        self.x -= 0.01

    if self.y <= ball.ymin_box + self.r: 
        self.vy = - self.vy
        self.y += 0.01
    elif self.y >= ball.ymax_box - self.r:
        self.vy = - self.vy
        self.y -= 0.01

    # for collisions
    for num2 in range(self.particle_num + 1, ball.Nparticles):
      # nx = (x2 - x1)
      nx = ball.particles[num2].x - self.x
      # ny = (y2 - y1)
      ny = ball.particles[num2].y - self.y
      # distance between center particles
      d = sqrt(nx**2 + ny**2)

      if d <= self.r + ball.particles[num2].r:

          # a is speed of particle 1, b the speed of particle 2
          ax = self.vx
          ay = self.vy
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
          an = ax*unx + ay*uny
          bn = bx*unx + by*uny

          at = ax*utx + ay*uty
          bt = bx*utx + by*uty

          # magnitude vectors after collision (f --> final)
          an_f = (an*(m1 - m2) + 2 * m2 * bn) / (m1 + m2)
          bn_f = (bn*(m2 - m1) + 2 * m1 * an) / (m1 + m2)
          at_f = at
          bt_f = bt

          # for particle 1
          self.vx = an_f * unx + at_f * utx
          self.vy = an_f * uny + at_f * uty

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
    # om er voor te zorgen dat bepaalde dingen in de update functie maar 1 keer gerunt worden per animate functie
    ball.counter = 0

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