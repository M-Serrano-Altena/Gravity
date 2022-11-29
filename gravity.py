# Marc Serrano Altena
# 29-11-2022
# this program simulates the gravitational pull of celestial body's

from matplotlib import pyplot as plt
from matplotlib import animation
from math import *
import random

class ball:
  dt = 0.01
  counter = 0

  # col stands for collision
  col_count_x = 0
  col_count_y = 0
  vx_col = 0
  vy_col = 0

  particle_num = 0
  xmin_box = 0
  xmax_box = 10
  ymin_box = 0
  ymax_box = 10

  particles = []
  Nparticles = 1
  circles = []

  fig, ax = plt.subplots()
  ax.set_xlim([xmin_box, xmax_box])
  ax.set_ylim([ymin_box, ymax_box])
  ax.set_aspect('equal')


  def __init__(self):
    self.x = random.randint(ball.xmin_box + 1, ball.xmax_box - 1) 
    self.y = random.randint(ball.ymin_box + 1, ball.ymax_box - 1)
    
    self.v = random.uniform(0.5, 1.5)
    self.r = random.uniform(0.15, 0.5)
    self.m = pi * self.r**2
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


  def update(self):
    ball.counter += 1
    collisionx = False
    collisiony = False

    # update atributes
    self.x += self.vx * ball.dt
    self.y += self.vy * ball.dt

    if self.vx != 0:
      self.vx += self.ax * ball.dt

    if self.vy != 0:
      self.vy += self.ay * ball.dt

    # orbits center of the box
    d_center_x = self.x - (ball.xmax_box - ball.xmin_box)/2
    d_center_y = self.y - (ball.ymax_box - ball.ymin_box)/2
    d_center = sqrt(d_center_x**2 + d_center_y**2)

    # unit vector
    self.uax = -d_center_x / d_center
    self.uay = -d_center_y / d_center

    if d_center > 1:
      self.a = 1/d_center**2

    self.ax = self.a * self.uax
    self.ay = self.a * self.uay

    # when hitting a wall; ± radius so that the bounce happens at the surface of the particle
    if self.x <= ball.xmin_box + self.r: 
        self.vx = - self.vx
        self.x += 0.01
        collisionx = True
    elif self.x >= ball.xmax_box - self.r:
        self.vx = - self.vx
        self.x -= 0.01
        collisionx = True

    if self.y <= ball.ymin_box + self.r: 
        self.vy = - self.vy
        self.y += 0.01
        collisiony = True
    elif self.y >= ball.ymax_box - self.r:
        self.vy = - self.vy
        self.y -= 0.01
        collisiony = True

    # collision with walls is not fully inelastic
    if collisionx:
      collisionx = False

      # so that a ball that should be at rest is actually at rest
      if ball.counter >= ball.col_count_x + 10 or ball.counter == ball.col_count_x:

        if round(self.vx, 2) == round(ball.vx_col, 2) and ball.counter >= ball.col_count_x + 4:
          self.vx = 0
          ball.counter = 0

        ball.col_count_x = ball.counter
        ball.vx_col = self.vx

      if self.vx > 0:
        self.vx -= 0.25

      elif self.vx < 0:
        self.vx += 0.25


    if collisiony:
      collisiony = False
      print(self.vy, ball.vy_col)

      # so that a ball that should be at rest is actually at rest
      if ball.counter >= ball.col_count_y + 50 or ball.counter == ball.col_count_y:

        if round(self.vy, 2) == round(ball.vy_col, 2) and ball.counter >= ball.col_count_y + 20:
          self.vy = 0
          ball.counter = 0

        ball.col_count_y = ball.counter
        ball.vy_col = self.vy

      if self.vy > 0:
        self.vy -= 0.25

      elif self.vx < 0:
        self.vy += 0.25

    # for collisions
    for num2 in range(self.particle_num + 1, ball.Nparticles):
      # nx = (x2 - x1)
      nx = ball.particles[num2].x - self.x
      # ny = (y2 - y1)
      ny = ball.particles[num2].y - self.y
      # distance between center particles
      d = sqrt(nx**2 + ny**2)

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


    # gravitational pull
    # self.ay += 0.1
          

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