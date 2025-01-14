import matplotlib.pyplot as plt
from scipy.constants import physical_constants
import numpy as np

def liikeyhtalo(x, v, k, m, c, L0):
    """
    Yksiulotteisen atomihilan liikeyhtälö.

    Parameters
    ----------
    x : Numpy array
        atomien paikat
    v : Numpy array
        atomien nopeudet
    k : float
        atomien välinen jousivakio
    m : float
        atomien massa
    c : float
        vaimennuskerroin
    L0 : float
        tasapainopituus

    Returns
    -------
    a : Numpy array
        atomien kiihtyvyys
    
    """
    
    a = np.zeros_like(x)
    for i in range(1, len(x) - 1):
        a[i] = k/m*(x[i - 1] - 2*x[i] + x[i + 1]) + c/m*(v[i - 1] - 2*v[i] + v[i + 1])
    
    a[0] = k/m*(x[1] - x[0] - L0)
    a[-1] = -k/m*(x[-1] - x[-2] - L0)
    
    return a

# Asetetaan tunnetut suureet muuttujiin
rho = 8570                  # kuparin tiheys [kg/m^3]
v_target = 5000             # haluttu äänen nopeus [m/s]
E = v_target**2 * rho       # Youngin moduli [Pa]
m = 63 * physical_constants['atomic unit of mass'][0]  # kupariatomin massa [kg]
N = 40                      # atomien määrä hilassa
L0 = (m / rho)**(1/3)       # kuparin yksikkökopin sivun pituus
k = E * L0                  # simulaatiossa käytettävä jousivakio kuparille
c = 0.01 * 2 * np.sqrt(m * k)  # vaimennuskerroin

# Alustetaan taulukot atomien paikkojen, nopeuksien ja 
# kiihtyvyyksien kirjanpitoa varten
x = L0 * np.arange(0, N)
v = np.zeros(N)
a = np.zeros(N)

# Alustetaan muut simulaatiossa tarvittavat suureet
T = 2 * np.pi / np.sqrt(k / m)  # Harmonisen oskillaattorin jaksonaika
dt = 0.01 * T                   # Aika-askel, joka on pieni verrattuna jaksonaikaan
dx = 2.5 * L0                   # Increase the initial displacement to ensure noticeable movement
x[0] = dx                       # Siirretään ketjun vasemmanpuoleista atomia dx verran
t = 0                           # Alustetaan aikalaskuri  
frame = 0                       # Alustetaan kuvien laskuri
max_images = 3                  # Maksimimäärä kuvia

# Tallennetaan äänen nopeus tiedostoon
sound_speed_data = []

# Päivitetään atomien paikkoja Eulerin menetelmän mukaisesti, 
# kunnes viimeinen atomi on liikahtanut 1% verran tasapainopaikaltaan
# eli kun "aalto" saavuttaa ketjun oikean reunan

while np.abs(x[-1] - L0*(N - 1)) < 0.01*dx:
    # Lisätään simulaation sisäistä aikaa
    t += dt
    
    # Debug print to check the values of t, x, and v
    print(f't = {t}, x = {x}, v = {v}')

    # Lasketaan liikeyhtälön avulla atomien kiihtyvyydet
    a = liikeyhtalo(x, v, k, m, c, L0)
    
    # Päivitetään kiihtyvyyden avulla uudet nopeudet, ja niiden avulla
    # uudet paikat atomeille
    for i in range(N):
        v[i] += dt * a[i]
        x[i] += dt * v[i]

    # Tallennetaan kuva joka 100. aika-askeleella
    if int(t/dt) % 100 == 0:
        plt.figure(figsize=(10, 6))
        plt.plot(x, np.zeros(N), 'bo')
        sound_speed = round((N - 1) * L0 / t, 2)
        plt.title(f'Atomien paikat ajanhetkellä, Sound Speed: {sound_speed} m/s')
        plt.xlabel('Paikka [m]')
        plt.ylabel('Atomit')
        plt.xlim(min(x) - 0.001 * L0, max(x) + 0.001 * L0)  # Further adjust x-axis limits to zoom in
        plt.ylim(-0.001, 0.001)  # Further adjust y-axis limits to zoom in
        plt.grid(True)
        plt.savefig(f'frame_{frame:04d}.png')
        plt.close()
        
        # Save sound speed to list
        sound_speed_data.append((t*1000, sound_speed))
        
        # Plot and save acceleration
        plt.figure(figsize=(10, 6))
        plt.plot(np.arange(N), a, 'ro-')
        plt.title(f'Atomien kiihtyvyys ajanhetkellä')
        plt.xlabel('Atomi')
        plt.ylabel('Kiihtyvyys [m/s^2]')
        plt.grid(True)
        plt.savefig(f'acceleration_{frame:04d}.png')
        plt.close()
        
        frame += 1

# Calculate the sound speed at the end of the simulation
sound_speed = round((N - 1) * L0 / t, 2)
print('Kuparin äänen nopeudeksi saatiin: ', sound_speed, '[m/s]')

# Display sound speed data in the terminal
print('Sound Speed (m/s)')
for time, speed in sound_speed_data:
    print(f'> {speed} <')

# Plot the sound speed over time
times, speeds = zip(*sound_speed_data)
plt.figure()
plt.plot(times, speeds, 'b-')
plt.title('Sound Speed Over Time')
plt.ylabel('Sound Speed (m/s)')
plt.xlabel('Time')
plt.grid(True)
plt.savefig('sound_speed_over_time.png')
#plt.show()

