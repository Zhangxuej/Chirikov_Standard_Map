import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

def chirikov_standard_map(K, x, p):
    '''
    input:  float, float
    output: float, float
    
    Returns the next iteration of the chrikov standard map. It 
    takes in initial x and p values and returns the subsequent
    x and p values given by the standard map:
    
    p_new = p + Ksin(x)
    x_new = x + p_new
    '''
    
    p_new = p + K*np.sin(x)
    x_new = x + p_new
    
    return (x_new, p_new)


def generate_all_Chirikov_maps(N):
    '''
    Generates the Chirikov standard map for a set number of points per orbit of N
    for K values running from 0 to 10 in increments of 0.02
    '''
    
    pi_const = 2*np.pi

    n = 0
    for K in np.arange(0, 10, 0.02):

        # Setting up figure
        fig, ax = plt.subplots(figsize=(10.0,10.0))

        title = "Chirikov Standard Map for K = {0}".format(str(round(K,3)))
        plt.title(title, fontsize=30)

        # Running over all x values
        for x_i in np.arange(0.0, 2*np.pi, 0.4):

            # Set up the arrays for x and p for N points
            x = np.zeros(N)
            p = np.zeros(N)

            # Setting initial value for x
            x[0] = x_i

            # Running over all initial p values from 0 to 2pi
            for p_i in np.arange(0.0, 2*np.pi, 0.4):

                # Setting initial value for p
                p[0] = p_i

                # Generating the next n-1 points
                for i in range(1,N):

                    # Apply Chirikov standard map function
                    x[i], p[i] = chirikov_standard_map(K, x[i-1], p[i-1])

                # Return values within map of 2pi x 2pi
                x = x % (2*np.pi)
                p = p % (2*np.pi)

                # Plot the result
                plt.scatter(x, p, s=0.0001, color='black')

        plt.xlabel("x/2pi", fontsize=30)
        plt.ylabel("p/2pi", fontsize=30)

        plt.xlim((0, 2*np.pi))
        plt.ylim((0, 2*np.pi))

        filename = "Chirikov_K_({0})".format(str(n))
        print(filename)
        plt.savefig(filename)

        n+=1


def generate_Chirikov_orbit(N, K, x0, p0):
    '''
    Generatures the Chirikov map for an individual (x0,p0) initial point
    '''
    
    # Set up the arrays for x and p for N points
    x = np.zeros(N)
    p = np.zeros(N)

    # Setting initial value for x
    x[0] = x0

    # Setting initial value for p
    p[0] = p0

    # Generating the next n-1 points
    for i in range(1,N):

        # Apply Chirikov standard map function
        x[i], p[i] = chirikov_standard_map(K, x[i-1], p[i-1])

    # Return values within map of 2pi x 2pi
    x = x % (2*np.pi)
    p = p % (2*np.pi)
    
    return x, p

 

# Number of points per standard map
N = 10000




#################################################################################
# Generating Chirikov Maps for Varying K for select (x,p) orbits
#################################################################################

'''
# Set up "kick" strength K
for K in np.arange(0.0, 5.0, 0.1):

    # Setting up figure
    fig, ax = plt.subplots(figsize=(20.0,20.0))

    title = "Chirikov Standard Map for K = {0}".format(str(K))
    plt.title(title, fontsize=30)


    # We generate the background image of the standard map upon which we will
    # plot our orbits over

    # Running over all x values
    for x_i in np.arange(0.0, 2*np.pi, 0.5):
        
        # Set up the arrays for x and p for N points
        x = np.zeros(N)
        p = np.zeros(N)
        
        # Setting initial value for x
        x[0] = x_i
        
        # Running over all initial p values from 0 to 2pi
        for p_i in np.arange(0.0, 2*np.pi, 0.5):
            
            # Setting initial value for p
            p[0] = p_i
            
            # Generating the next n-1 points
            for i in range(1,N):
                
                # Apply Chirikov standard map function
                x[i], p[i] = chirikov_standard_map(K, x[i-1], p[i-1])
            
            # Return values within map of 2pi x 2pi
            x = x % (2*np.pi)
            p = p % (2*np.pi)
            
            # Plot the result
            plt.scatter(x, p, s=0.00001, color='black')


    # Generate the orbits for initial (x,p) pairs:
    # (0.5,0.5), (0.0, 3.0), and (2pi/3, pi)
    orbit_x_1, orbit_p_1 = generate_Chirikov_orbit(N, K, 0.5, 0.5)
    orbit_x_2, orbit_p_2 = generate_Chirikov_orbit(N, K, 0.0, 3.0)
    orbit_x_3, orbit_p_3 = generate_Chirikov_orbit(N, K, (2/3)*np.pi, np.pi)

    # Plot these orbits as scatter graphs with different colours
    plt.scatter(orbit_x_1, orbit_p_1, s=0.1, color='blue')
    plt.scatter(orbit_x_2, orbit_p_2, s=0.1, color='red')
    plt.scatter(orbit_x_3, orbit_p_3, s=0.1, color='green')
    plt.xlabel("x mod 2pi", fontsize=30)
    plt.ylabel("p mod 2pi", fontsize=30)


    # Limit our x and y ranges to 0 to 2pi
    plt.xlim((0, 2*np.pi))
    plt.ylim((0, 2*np.pi))

    # Give each map its own image file and filename
    filename = "Chirikov_Orbit_K_({0}).png".format(str(K))
    print(filename)
    plt.savefig(filename)

'''

#################################################################################
# Generating Orbit as a Function of Time for Chirikov Map with set K
#################################################################################


# We want to compare the time evolution of the orbits for two particular values of
# K: K_1 exhibits normal periodic motion whereas K_2 exhibits chaotic motion
# K_3 showcases periodic, quasi-periodic, and chaotic motion altogether which makes
# it a particularly interesting sample to draw plots over time

# For each map we will take a look at the orbits with initial points (0.5,0.5), 
# (0.0, 3.0), and (2pi/3, pi). We will take a look at these individually
#initial_point = (0.5, 0.5)
#initial_point = (0.0, 3.0)
initial_point = ((2.0*np.pi)/3, np.pi)

# Define our three maps we want to plot our orbits for
K_1 = 0.02
K_2 = 3.0
K_3 = 0.971635


# We will store the points into the an array for each K value
time_orbit = np.zeros((N,2))
time_orbit[0] = initial_point

# Iterate over each initial point
for n in range(1, 200):

    # Setting up figure
    fig, ax = plt.subplots(figsize=(10.0,10.0))

    # Calculate the new instance of the orbit at t = n 
    time_orbit[n] = chirikov_standard_map(K_3, time_orbit[n-1][0], time_orbit[n-1][1])

    # Setting up title of the plot
    title = "Orbit of Initial Point ({0}, {1}) as a Function of Time, K = {2}".format(str(round(initial_point[0],4)), str(round(initial_point[1], 4)), str(K_3)  )
    plt.title(title, fontsize=20)

    # Running over all x values
    for x_i in np.arange(0.0, 2*np.pi, 0.5):

        # Set up the arrays for x and p for N points
        x = np.zeros(N)
        p = np.zeros(N)

        # Setting initial value for x
        x[0] = x_i

        # Running over all initial p values from 0 to 2pi
        for p_i in np.arange(0.0, 2*np.pi, 0.5):

            # Setting initial value for p
            p[0] = p_i

            # Generating the next n-1 points
            for i in range(1,N):

                # Apply Chirikov standard map function
                x[i], p[i] = chirikov_standard_map(K_3, x[i-1], p[i-1])

            # Return values within map of 2pi x 2pi
            x = x % (2*np.pi)
            p = p % (2*np.pi)

            # Plot the result
            plt.scatter(x, p, s=0.00001, color='black')

    # We make sure to apply the modulus of 2pi before plotting the points
    plt.scatter(time_orbit[n][0] % (2*np.pi), time_orbit[n][1] % (2*np.pi), s=10.0, color='blue')

    plt.xlabel("x mod(2π)", fontsize=30)
    plt.ylabel("p mod(2π)", fontsize=30)

    plt.xlim((0, 2*np.pi))
    plt.ylim((0, 2*np.pi))

    filename = "{0}".format(str(n))
    print(filename)
    plt.savefig(filename)





