from pylab import *
from scipy.integrate import *

def rot_saddle_trap_ode_x_y(t,x,q=0.5):
    # x = x , vx , y , vy
    dx , dy = x[1] , x[3]
    dvx = -2*q*(x[0]*cos(2*t) + x[2]*sin(2*t))
    dvy = 2*q*(x[2]*cos(2*t) - x[0]*sin(2*t))
    return array([dx,dvx,dy,dvy])

def nu2q(nu):
    g , h , R = 9810 , 30 , 54 # default parameters for the trap
    p=(g*h)/(R*R)
    omega = 2*pi*nu
    return p/(omega**2)

def draw_simulation_trajectory(x_sol,y_sol):
    fig, ax = plt.subplots(dpi=150)
    ax.plot(x_sol,y_sol)
    ax.grid(),title("Trajectory"),xlabel(r"$x$ [mm]"),ylabel(r"$y$ [mm]")
    ax.set_aspect('equal')
    show()

def initial_cond():
    print("Welcome \U0001F600 \n")
    print('Enter the initial values for x and y')
    x_ini = float(input('\t x = '))
    y_ini = float(input('\t y = '))
    return x_ini , y_ini

def run_simulation():
    x_ini,y_ini = initial_cond()
    ini_p = [x_ini,y_ini,0.0,0.0] #####
    t = (0,500)
    my_q = nu2q(4)
    xy_sol = solve_ivp(rot_saddle_trap_ode_x_y,t,ini_p,args=(my_q,),method='LSODA',dense_output=True,rtol=1e-8,atol=1e-8)
    if xy_sol.success == True:
        print("I have the solution for the entered values")
    else:
        print("I don't have the solution for the entered values")
    t_cont = linspace(0,xy_sol.t[-1],5000)
    x_sol , y_sol = xy_sol.sol(t_cont)[0,:] , xy_sol.sol(t_cont)[2,:]
    draw_simulation_trajectory(x_sol , y_sol)

#run_simulation()
