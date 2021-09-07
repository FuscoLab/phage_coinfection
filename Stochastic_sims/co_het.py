import random
import numpy as np
import time
import math
import pickle
import os
import copy
import sys
# import matplotlib.pyplot as plt

time_start = time.time()

time_step=0
max_time = 120000
initial_virus = 0
initial_bacteria = int(sys.argv[2])
initial_virus_2 = int(initial_bacteria*2)
A = 1
y_2 = float(sys.argv[3])
y_1 = float(sys.argv[3])
lyse_time_2 = int(sys.argv[5])
lyse_time_1 = int(sys.argv[5])
alpha = float(sys.argv[4])
beta = float(sys.argv[4])
delta = 0.1

repeat_main=int(sys.argv[6])
steady_time=500 #time to (stochastically) define the steady state
t_back=30





class Bacteria:
    def __init__(self, location):
        self.location = location
        self.infection_time = 0
        self.alive = 1
        self.virus_1 =0
        self.virus_2= 0
        self.primary = 0
    def lyse(self,output_1,output_2,time_step,count_bacteria_alive_bin,j,count_1_infect,count_2_infect):
        if self.primary==2:
            if self.infection_time+lyse_time_2==time_step:
                burst_size=np.random.poisson(y_2)
                burst_v1=np.random.binomial(burst_size,self.virus_1/(self.virus_1+self.virus_2))
                output_1+=burst_v1
                output_2+=burst_size-burst_v1
                count_2_infect-=1
                self.infection_time = 0
                self.alive = 1
                self.virus_1 =0
                self.virus_2= 0
                self.primary = 0
        elif self.primary==1:
            if self.infection_time+lyse_time_1==time_step:
                burst_size=np.random.poisson(y_1)
                burst_v1=np.random.binomial(burst_size,self.virus_1/(self.virus_1+self.virus_2))
                output_1+=burst_v1
                output_2+=burst_size-burst_v1
                count_1_infect-=1
                self.infection_time = 0
                self.alive = 1
                self.virus_1 =0
                self.virus_2= 0
                self.primary = 0
        else:
            count_bacteria_alive_bin +=1
        return count_bacteria_alive_bin,output_1,output_2,count_1_infect,count_2_infect

def infection(bacteria_set,virus_1,virus_2,time_step,count_bacteria_infected,count_bacteria_alive,count_1_infect,count_2_infect):
    number_alive=initial_bacteria
    count_bacteria_infected=0
    for position in range(int(number_alive)):
        if bacteria_set[position].primary!=0:
            count_bacteria_infected+=1
    Y= np.random.poisson(beta*virus_2*number_alive)
    X= np.random.poisson(alpha*virus_1*number_alive)
    X=np.clip(X,a_max=virus_1,a_min=0)
    Y=np.clip(Y,a_max=virus_2,a_min=0)
    indices_2 = np.random.randint(number_alive,size=(int(Y))) if number_alive != 0 else np.empty(0)
    indices_1 = np.random.randint(number_alive,size=(int(X))) if number_alive != 0 else np.empty(0)
    random_prob = np.random.random(len(np.intersect1d(indices_1,indices_2)))
    countX_intersect= np.array([])
    countY_intersect= np.array([])
    intersections = np.intersect1d(indices_1,indices_2)
    for i in range(len(intersections)):
        countX_intersect=np.append(countX_intersect, np.count_nonzero(indices_1==intersections[i]))
        countY_intersect=np.append(countY_intersect, np.count_nonzero(indices_2==intersections[i]))
    prob=countX_intersect/(countX_intersect+countY_intersect)
    virus_1-=len(indices_1)
    virus_2-=len(indices_2)
    indices_1_extra=np.array([])
    indices_2_extra=np.array([])
    for i in range(len(intersections)):
        if random_prob[i]>prob[i]:
            indices_1_extra=np.append(indices_1_extra, indices_1[(indices_1==intersections[i])])
            indices_1=indices_1[~(indices_1==intersections[i])]
        elif random_prob[i]<prob[i]:
            indices_2_extra=np.append(indices_2_extra, indices_2[(indices_2==intersections[i])])
            indices_2=indices_2[~(indices_2==intersections[i])]
        else:
            coinflip=np.random.binomial(1, 0.5)
            if coinflip==0:
                indices_1_extra=np.append(indices_1_extra, indices_1[(indices_1==intersections[i])])
                indices_1=indices_1[~(indices_1==intersections[i])]
            else:
                indices_2_extra=np.append(indices_2_extra, indices_2[(indices_2==intersections[i])])
                indices_2=indices_2[~(indices_2==intersections[i])]
    indices_1_extra=indices_1_extra.astype(int)
    indices_2_extra=indices_2_extra.astype(int)
    for position in indices_1:
        bacteria_set[position].virus_1+=1
        if bacteria_set[position].infection_time==0:
            bacteria_set[position].infection_time=time_step
            bacteria_set[position].primary=1
            count_bacteria_infected+=1
            count_1_infect+=1
    for position in indices_2:
        bacteria_set[position].virus_2+=1
        if bacteria_set[position].infection_time==0:
            bacteria_set[position].infection_time=time_step
            bacteria_set[position].primary=2
            count_bacteria_infected+=1
            count_2_infect+=1
    for position in indices_1_extra:
        bacteria_set[position].virus_1+=1
    for position in indices_2_extra:
        bacteria_set[position].virus_2+=1
    return bacteria_set,virus_1,virus_2,count_bacteria_infected,count_1_infect,count_2_infect

def viral_decay(virus_1_time,virus_2_time):
    virus_1_time-=np.random.poisson(virus_1_time*delta)
    virus_1_time=np.clip(virus_1_time,a_min=0,a_max=None)
    virus_2_time-=np.random.poisson(virus_2_time*delta)
    virus_2_time=np.clip(virus_2_time,a_min=0,a_max=None)
    return virus_1_time,virus_2_time

def viral_growth(bacteria_set):
    for position in range(initial_bacteria):
        if bacteria_set[position].primary==1:
            temp_v1=float((y_1/lyse_time_1)*bacteria_set[position].virus_1/(bacteria_set[position].virus_1+bacteria_set[position].virus_2))
            temp_v2=float((y_1/lyse_time_1)*bacteria_set[position].virus_2/(bacteria_set[position].virus_1+bacteria_set[position].virus_2))
            bacteria_set[position].virus_1+=temp_v1
            bacteria_set[position].virus_2+=temp_v2
        elif bacteria_set[position].primary==2:
            temp_v1=float((y_2/lyse_time_2)*bacteria_set[position].virus_1/(bacteria_set[position].virus_1+bacteria_set[position].virus_2))
            temp_v2=float((y_2/lyse_time_2)*bacteria_set[position].virus_2/(bacteria_set[position].virus_1+bacteria_set[position].virus_2))
            bacteria_set[position].virus_1+=temp_v1
            bacteria_set[position].virus_2+=temp_v2
    return bacteria_set

np.random.seed(10*int(sys.argv[1]))
heterozygosity_sum=np.zeros(max_time)
heterozygosity_nonzero_sum=np.zeros(max_time)

# growth_rates=np.zeros((repeat_main,1))
# intercept=np.zeros((repeat_main,1))
for z in range(repeat_main):
    time_step=0
    move_flag=0
    fix_tester=0
    move_counter = 0
    bacteria_set = [Bacteria(0) for j in range(initial_bacteria)]
    virus_1 = np.zeros(max_time)
    virus_1[0] = initial_virus
    virus_2 = np.zeros(max_time)
    virus_2[0] = initial_virus_2
    heterozygosity=np.zeros(max_time)
    heterozygosity_nonzero=np.zeros(max_time)
    heterozygosity[0]=2*virus_1[0]*virus_2[0]/(virus_1[0]+virus_2[0])
    heterozygosity_nonzero[0]=heterozygosity[0]>0
    x = np.arange(max_time)
    count_bacteria_alive = np.zeros(max_time,dtype=int)
    count_bacteria_infected = np.zeros(max_time,dtype=int)
    count_bacteria_alive[0] = initial_bacteria
    count_1_infect=0
    count_2_infect=0
    time_step+=1
    steady = steady_time+np.random.randint(0,600)
    # print(steady)

    #Generates the steady state starting population?
    while (time_step < steady):
        virus_1[time_step]=virus_1[time_step-1]
        virus_2[time_step]=virus_2[time_step-1]
        output_1=0
        output_2=0
        bacteria_set = viral_growth(bacteria_set)
        for j in range(initial_bacteria):
            if bacteria_set[j].alive == 1:
                count_bacteria_alive[time_step],output_1,output_2,count_1_infect,count_2_infect = bacteria_set[j].lyse(output_1,output_2,time_step,count_bacteria_alive[time_step],j,count_1_infect,count_2_infect)
        virus_1[time_step]+=output_1
        virus_2[time_step]+=output_2
        bacteria_set , virus_1[time_step] , virus_2[time_step],count_bacteria_infected[time_step],count_1_infect,count_2_infect = infection(bacteria_set,virus_1[time_step],virus_2[time_step],time_step,count_bacteria_infected[time_step],count_bacteria_alive[time_step],count_1_infect,count_2_infect)
        virus_1[time_step],virus_2[time_step]=viral_decay(virus_1[time_step],virus_2[time_step])
        heterozygosity[time_step]=2*virus_1[time_step]*virus_2[time_step]/((virus_1[time_step]+virus_2[time_step])**2)
        heterozygosity_nonzero[time_step]=heterozygosity[time_step]>0
        time_step +=1

    steady_size=np.mean(virus_2[(steady-100):steady])
    steady_size_error=np.std(virus_2[(steady-100):steady])/10
    virus_2_initiate=virus_2
    virus_1_initiate=virus_1
    count_1_infect_initiate=count_1_infect
    count_2_infect_initiate=count_2_infect
    count_bacteria_alive_initiate=count_bacteria_alive
    count_bacteria_infected_initiate=count_bacteria_infected
    bacteria_set_initiate=copy.deepcopy(bacteria_set)
    bacteria_set=copy.deepcopy(bacteria_set_initiate)
    virus_2=virus_2_initiate.copy()
    virus_1=virus_1_initiate.copy()
    count_bacteria_infected=count_bacteria_infected_initiate.copy()
    count_bacteria_alive=count_bacteria_alive_initiate.copy()
    count_1_infect=copy.copy(count_1_infect_initiate)
    count_2_infect=copy.copy(count_2_infect_initiate)
    time_step=steady
    while (time_step < max_time):
        virus_1[time_step]=virus_1[time_step-1]
        virus_2[time_step]=virus_2[time_step-1]
        if time_step==steady+1:
            mutant_fraction=0.5
            mutant_infected=int(mutant_fraction*count_bacteria_infected[time_step-1])
            virus_1[time_step]+=virus_2[time_step]/2
            virus_2[time_step]-=virus_2[time_step]/2
            for j in range(initial_bacteria):
                if bacteria_set[j].primary==2:
                    if mutant_infected!=0:
                        bacteria_set[j].primary=1
                        bacteria_set[j].virus_1=bacteria_set[j].virus_2
                        bacteria_set[j].virus_2=0
                        count_1_infect+=1
                        count_2_infect-=1
                        mutant_infected-=1
        output_1=0
        output_2=0
        bacteria_set = viral_growth(bacteria_set)
        for j in range(initial_bacteria):
            if bacteria_set[j].alive == 1:
                count_bacteria_alive[time_step],output_1,output_2,count_1_infect,count_2_infect = bacteria_set[j].lyse(output_1,output_2,time_step,count_bacteria_alive[time_step],j,count_1_infect,count_2_infect)
        virus_1[time_step]+=output_1
        virus_2[time_step]+=output_2
        bacteria_set , virus_1[time_step] , virus_2[time_step],count_bacteria_infected[time_step],count_1_infect,count_2_infect = infection(bacteria_set,virus_1[time_step],virus_2[time_step],time_step,count_bacteria_infected[time_step],count_bacteria_alive[time_step],count_1_infect,count_2_infect)
        virus_1[time_step],virus_2[time_step]=viral_decay(virus_1[time_step],virus_2[time_step])
        heterozygosity[time_step]=2*virus_1[time_step]*virus_2[time_step]/((virus_1[time_step]+virus_2[time_step])**2)
        heterozygosity_nonzero[time_step]=heterozygosity[time_step]>0
        time_step+=1
    print(repeat_main-z)
    heterozygosity_sum+=heterozygosity
    heterozygosity_nonzero_sum+=heterozygosity_nonzero

result=np.zeros((max_time, 2))
result[:,0]=heterozygosity_sum/repeat_main
result[:,1]=heterozygosity_nonzero_sum/repeat_main

np.savetxt('Het_co_'+str(y_1)+'burst'+str(y_2)+'_'+str(alpha)+'alpha'+str(beta)+'_'+str(lyse_time_1)+'tau'+str(lyse_time_2)+'_'+str(sys.argv[1])+'.txt',result)
