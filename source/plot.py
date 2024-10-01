import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

data1 = pd.read_csv("../data/s_1000.txt", sep='\s+', header=None)
data1 = pd.DataFrame(data1)

s1 = data1[0]

data2 = pd.read_csv("../data/s.txt", sep='\s+', header=None)
data2 = pd.DataFrame(data2)

s = data2[0]

data3 = pd.read_csv("../data/eigenvalues_1000.txt", sep='\s+', header=None)
data3 = pd.DataFrame(data3)

eigen1 = data3[0]

data4 = pd.read_csv("../data/eigenvalues.txt", sep='\s+', header=None)
data4 = pd.DataFrame(data4)

eigen = data4[0]

#########################################################################
#########################################################################
N = 700
beta = 100
sigma_diag = np.sqrt(1/(2*beta));

def semi_circle(x, a, b):
    bl = (1/(2*np.pi*b*b))*np.sqrt(4*a - x*x);
    return bl;

eigen1000 = []
for item in eigen1:
    if abs(item) < 1.5:
        eigen1000.append(item)

deltaE = []
for i in range(0, len(eigen1000) - 1):
    deltaE.append(abs(eigen1000[i] - eigen1000[i+1]))

s1000 = []
for j in range(0, len(deltaE)):
    s1000.append(deltaE[j]*semi_circle(eigen1000[j] ,N , sigma_diag))
    
print(len(s1000))


#########################################################################
#########################################################################

plt.subplot(1,2,1)

# Plotting histogram
hist, bins = np.histogram(eigen1000, bins=100, density=True)  # Density=True for normalized histogram



# Overlaying Gaussian function
x = np.linspace(-2, 2, 1000)
bin_centers = (bins[:-1] + bins[1:]) / 2

popt, pcov = curve_fit(semi_circle, bin_centers, hist)

# Normalizing the function
#integral = np.trapz(y, x)
#y_normalized = (y / integral)

#plt.plot(x, y, '--', label='Semi-cirlce law', color='red')
plt.hist(eigen1000, bins=100, edgecolor='black', density=True, alpha=0.7)
plt.plot(bin_centers, semi_circle(bin_centers, *popt), 'r-', label='Fitted Semi-circle law for eigenvalues')
plt.title('Histogram with Semi-circle law')
plt.xlabel('Eigenvalues')
plt.ylabel('Frequency')
plt.legend()

#########################################################################
#########################################################################

plt.subplot(1,2,2)

def wigner(x, a, b):
    bl = a*x*np.exp(-b*(np.pi/4)*(x)**2)
    return bl;

# Plotting histogram
hist, bins = np.histogram(s1, bins=100, density=True)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Overlaying Gaussian function
x = np.linspace(0, 40, 1000)
#y = (np.pi/2)*x*np.exp(-(np.pi/4)*x**2)
popt, pcov = curve_fit(wigner, bin_centers, hist)

# Normalizing the function
#integral = np.trapz(y, x)
#y_normalized = (y / integral)

#plt.plot(x, y_normalized, '--', label='Semi-cirlce law', color='red')

plt.hist(s1, bins=100, edgecolor='black', density=True, alpha=0.7)  # Density=True for normalized histogram
plt.plot(x, wigner(x, *popt), 'r-', label='Fitted Wigners formula for spacing probability (GOE)')

plt.title('Wigner formula')
plt.xlabel('Spacing (s)')
plt.ylabel('Frequency')
plt.legend()
plt.xlim(0, 40)

plt.show()