"""

@author: Alex Kerr
"""

import numpy as np

#constants
G = 1.

class StarSystem:
    """A collection of planets
    
    Args:
        pos (list-like): A list of lists of the positions of the planets, in AU's.
        vel (list-like): A list of lists of the initial velocities of the planets, in AU's/yr., indexed like pos.
        mass (list-like): A list of the masses of the planets in the system, in earth masses, indexed like pos.
        
    Keywords:
        force (str): String that selects the force interaction.
        dt (float): The time step size in integration.
        method(str): String that selects integration method."""
    
    def __init__(self, pos, mass, force="newton", dt=0.1, method="euler"):
        self.pos = np.array(pos)
        self.mass = np.array(mass)
        self._configure_interactions()
        
    def move(self):
        """Change the positions and velocities of the planets along a time step."""
        pos_step, vel_step = method_dict[self.method]
        self.pos += pos_step
        self.vel += vel_step
        
    def define_force_routine(self):
        """A factory function that returns the function that calculates the gradient on the star system."""
        
        i, j = self.iacts[:,0], self.iacts[:,1]
        posij = self.pos[i] - self.pos[j]
        rij = np.linalg.norm(posij, axis=1)
        
        if self.force == "newton":
            
            def grad_newton():
                grad = np.zeros((len(self.pos), 3))
                newton = G*(self.mass[i]*self.mass[j]/(rij**3))[:,None]*posij
                np.add.at(grad, i, newton)
                np.add.ad(grad, j, -newton)
                return grad
                
            self.calculate_force = grad_newton
            
        else:
            raise ValueError("Invalid force")
        
    def _configure_interactions(self):
        """Assign the interactions in the star system."""
        self.iact = np.array([[i,j] for i in range(len(self.pos)) for j in range(len(self.pos)) if  i < j])
        
        
def euler(system):
    return system.vel*system.dt, system.calculate_force()*system.dt/system.mass
    
method_dict = {"euler":euler}
        
tattooine = StarSystem([[0.,0.],[1.,1.]], [1., 0.9])
print(tattooine.iact) 