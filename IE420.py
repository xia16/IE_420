import numpy as np
import time
start_time = time.time()

#initialization

Option = raw_input('Enter C for Call and P for puts\n')
while (Option != 'C') and (Option != 'P'):
    print("Input Not Accepted")
    Option = raw_input('Enter C for Call and P for puts\n')
print("Option Input Accepted")

excercise = raw_input("Enter A for American and E for European\n")
while ((excercise != 'A') and (excercise != 'E')):
    print("Input Not Accepted")
    excercise = raw_input("Enter A for American and E for European\n")
print("EURO US Input Accepted\n")

K = int(raw_input("Enter K for strike price\n"))
T = int(raw_input("Enter T for time to maturity\n"))
S_0 = int(raw_input("Enter Initial stock price\n"))
V = float(raw_input("Enter Volatility\n"))
r = float(raw_input("Enter continous compounding risk free interest rate as %\n"))
q = float(raw_input("Enter continous dividend yield as %\n"))
N = int(raw_input("Enter Number of time steps\n"))
delta = int(T)/int(N)
print(delta)
u = np.exp(V*(delta)**(1/2))
d = np.exp(-V*(delta)**(1/2))
p_star = ((np.exp(delta*(r-q)-d)/(u-d)))

#Array for terminal stock price same for put and calls , A or E
terminal_SP = np.zeros((N+1,N+1))
terminal_SP[0][0] = S_0
for i in range(1,N+1):
    terminal_SP[i][0] = terminal_SP[i-1,0]*u
    for j in range(1,i+1):
        terminal_SP[i,j] = terminal_SP[i-1,j-1]*d
print(terminal_SP)

if excercise == 'E':
    if Option == 'P':
        temp1 = -1
        optval = np.zeros((N + 1, N + 1))
        for i in range (N+1):
            optval[N][i] = max(0,temp1*(terminal_SP[N][i]-K))
        for i in range(N - 1, -1, -1):
            for j in range(N - 1, -1, -1):
                optval[i][j] = (p_star * optval[i - 1][j] + (1 - p_star) * optval[i - 1][j + 1])
    if Option == 'C':
        temp1 = 1
        optval = np.zeros((N + 1, N + 1))
        for i in range(N + 1):
            optval[N][i] = max(0, temp1 * (terminal_SP[N][i] - K))
        for i in range(N-1,-1,-1):
            for j in range(N-1,-1,-1):
                optval[i][j] = (p_star*optval[i-1][j]+(1-p_star)*optval[i-1][j+1])
if excercise == 'A':
    if Option == 'P':
        temp1 = -1
        optval = np.zeros((N + 1, N + 1))
        for i in range (N+1):
            optval[N][i] = max(0,temp1*(terminal_SP[N][i]-K))
        for i in range(N-1,-1,-1):
            for j in range(N-1,-1,-1):
                optval[i][j] = max(temp1*(terminal_SP[i][j]-K),p_star*optval[i-1][j]+(1-p_star)*optval[i-1][j+1])
    if Option == 'C':
        temp1 = 1
        optval = np.zeros((N + 1, N + 1))
        for i in range(N + 1):
            optval[N][i] = max(0, temp1 * (terminal_SP[N][i] - K))
        for i in range(N-1,-1,-1):
            for j in range(N-1,-1,-1):
                optval[i][j] = max(temp1*(terminal_SP[i][j]-K),p_star*optval[i-1][j]+(1-p_star)*optval[i-1][j+1])
print(optval[0][0])
print("--- %s seconds ---" % (time.time() - start_time))