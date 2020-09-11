# -*- coding: utf-8 -*-
import random
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd
import numpy as np
from scipy import stats


'''
Global variables Definition
'''
index = 0 #The index of each case
casesList = [] #Save the details of each case
last_day = 132 #Simulation period，the first day is Jan 1, 2020
sim_num = 10 #The number of simulation times of each scenario
asymptomatic_inf = 0.1 #Relative Infectiousness of asymptomatic cases


#get tp
def get_tracing_rate(parent_confirmed_time):
    r = 0
    
    if parent_confirmed_time < 54:#3月16日
        r = 0.31
    elif parent_confirmed_time < 76:#4月7日
        r = 0.13
    elif parent_confirmed_time < 132:#4月7日
        r = 0.5
    else:
        r = 0.5
    
    xk = (0,1)
    pk = (1-r, r)
    
    return int(stats.rv_discrete(values=(xk, pk)).rvs())


#get ap
def get_awareness_rate(parent_confirmed_time):
    r = 0
    
    if parent_confirmed_time < 54:#Mar 16
        r = 0.54
    elif parent_confirmed_time < 76:#Apr 7
        r = 0.25
    elif parent_confirmed_time < 132:#Jun 1
        r = 0.4
    else:
        r = 0.4
        
    xk = (0,1)
    pk = (1-r, r)
    
    return int(stats.rv_discrete(values=(xk, pk)).rvs())



#get proportion of asymptomatic cases
def get_is_asymptomatic():
    r = 0.15
    xk = (0,1)
    pk = (1-r, r)
    
    return int(stats.rv_discrete(values=(xk, pk)).rvs())


#get secondary number of each caes according to the infect time
def get_secondary_num(infect_time, source,is_asymptomatic):

    R = 0
    
    if source == 'Local':
        if infect_time < 54:#3月16日
            R = 2.5
        elif infect_time < 76:#4月7日
            R = 1.51
        elif infect_time < 132:#6/1
            R = 1.29
        else:
            R = 2.5
            
    else:
        if infect_time < 54:#3月16日
            R = 2.5
            
        elif infect_time < 63:#3月25
            R = 1.51
        
        else:
            R = 0
    
    if is_asymptomatic == 1:
        R = R *asymptomatic_inf
    
    if R <= 0:
        return 0
    else:
        
        p = 0.48

        n=(p*R)/(1-p)
        #return min(10,max(0,int(st.nbinom.rvs(n=n,p=p,size=1,loc=0))))
        return max(0,int(st.nbinom.rvs(n=n,p=p,size=1,loc=0)))
    
#infection to symptom onset delay for imported cases
def import_infect2sym():
    return min(7,max(0,int(round(st.skewnorm.rvs(a=2,loc=-2,scale=4)))))

#onset to isolation delay for imported cases
def distr_logistic():
    #return min(23,max(0,int(round(st.t.rvs(df=1.38, loc=0.91, scale=0.96)))))
    #return round(st.skewnorm.rvs(a=18855325.0091, loc=-0.0000, scale=3.1015))
    #return int(st.wald.rvs(loc=-0.4964, scale=2.2623))
    #return st.foldcauchy.rvs(c=0.0002, loc=0, scale=0.9999)
    
    r = random.random()*100
    A1 = 0.37899
    A2 = 0.02263
    x0 = 1.41316
    p = 3.1005
    
    temp = 0
    for x in range(11):
        y = 100*(A2 + (A1-A2)/(1+(x/x0)**p))   
        temp += y
        if r <= temp:
            return x
    
    


#get infect time of new cases 
def get_infected_time(infect2sym_time,size):
    return st.skewnorm.rvs(a=10, scale=3,loc=infect2sym_time,size=size)
    #return st.skewnorm.rvs(a=2, scale=3,loc=infect2sym_time,size=size)
    
#get incubation period
def incubation_period():  
    return int(st.weibull_min.rvs(c=1.88,loc=0,scale=7.97))
    #return min(20,int(st.weibull_min.rvs(c=1.88,loc=0,scale=7.97)))
    #return int(distr_exponweib(2.322737,6.492272,1,loc=0,scale=5))
    #return max(0,int(round(st.exponweib.rvs(a=2.322737, c=6.492272,loc=0,scale=5))))

#serial interval
def serial_interval():
    return round(st.skewnorm.rvs(a=1.1122, loc=0.3588, scale=5.8032))
    
#tracing delay
def tracing_isolation_time():
    return min(4,round(st.skewnorm.rvs(a=5818014.6979, loc=0, scale=2.3836)))
    
#isolation to confirmation delay 
def confirmed_delay():
    return min(6,max(0,round(st.lognorm.rvs(s=0.6110, loc=-0.5353, scale=1.6279))))

#onset to isolation delay for class3a
def sym2isolation_delay_3a():
    #return min(9,int(round(st.skewnorm.rvs(a=15958514.68, scale=7,loc=0))))
    return st.exponweib.rvs(a=0.1943, c=2.1613, loc=0, scale=15.3507)

#onset to isolation delay for class3b
def sym2isolation_delay_3b():
    return min(8,int(round(st.foldcauchy.rvs(c=1.3474, loc=0, scale=1.7909))))

#class4 发病到隔离
def sym2isolation_delay_4(): 
    return min(17,int(round(st.lognorm.rvs(s=0.5792, loc=0.6142, scale=6.0806))))

#recovered time, here we didn't take into account
def recovered_time():
    #return int(st.mielke.rvs(k=2.0670, s=5.2403, loc=-0.4848, scale=15.7656))
    return 0
    



#create new case according to parent case (infector)
def create_new_case(parent_case):
    
    #如果输入病例别被隔离
    if parent_case[1] == 0 and parent_case[8] == '1':
        return
    
    
    if parent_case[2] >= last_day: # or parent_case[2]<-10
        return
    
    if parent_case[3] >= last_day:
        return
    
    
    parent_infected_time = parent_case[2]
    parent_sym_time = parent_case[3]
    parent_isolation_time = parent_case[4]
    parent_confirmed_time = parent_case[5]
    parent_source = parent_case[7]
    parent_is_asymptomatic = parent_case[9]

    new_cases_num = get_secondary_num(parent_infected_time,parent_source,parent_is_asymptomatic)
    
    for i in range(new_cases_num):
        
        for j in range(100):
            
            sym_time = parent_sym_time + serial_interval()
            infected_time = sym_time - incubation_period()
            if sym_time > parent_infected_time and infected_time > parent_infected_time:
                break
        '''   
        for k in range(100):
            infected_time = sym_time- incubation_period()#感染时间
            if infected_time > parent_infected_time:
                break
        '''
        if sym_time < parent_infected_time:
            sym_time = parent_infected_time
        if infected_time < parent_infected_time:
            infected_time = parent_infected_time
        
        if parent_is_asymptomatic == 0 and infected_time >= parent_isolation_time:
            break
        
        
        global index
        index += 1 
        parent_index = parent_case[0]  
        
        
        isolation_time = 0 
        confirmed_time = 0 
        recover_time = 0 
        
        source = 'Local'
        class_type = ''
        

        is_asymptomatic = get_is_asymptomatic()
        

        isolation_time = parent_confirmed_time + tracing_isolation_time()
        if is_asymptomatic == 1:
            isolation_time = -100
            class_type = '-100'
            confirmed_time = -100
        
        elif get_tracing_rate(parent_confirmed_time) == 1 and isolation_time-sym_time < 10:
            
            if sym_time < isolation_time:
                class_type = '2'
                confirmed_time = isolation_time + confirmed_delay()
                
            else:
                class_type = '1'
                confirmed_time = sym_time + confirmed_delay() 
        

        else:
            
            if get_awareness_rate(parent_confirmed_time) == 1:
                if sym_time < parent_case[5]:
                    class_type = '3a'
                    isolation_time = sym_time + sym2isolation_delay_3a() 
                    
                else:
                    class_type = '3b'
                    isolation_time = sym_time + sym2isolation_delay_3b()
                    
                    
            else:
                class_type = '4'
                isolation_time = sym_time + sym2isolation_delay_4()
            
            
            confirmed_time = isolation_time+ confirmed_delay()
        
        
        recover_time = confirmed_time + recovered_time()
        
        case = [index, parent_index, infected_time, sym_time, isolation_time, confirmed_time, recover_time, source, class_type,is_asymptomatic]
        
        if case[2]>last_day:
            return
        
        casesList.append(case)
        create_new_case(case)
            
                
      

#create imported cases
def imported_case(today):
    global index
    index += 1
    parent_index = 0 
    confirmed_time = today 
    isolation_time = 0
    sym_time = 0
    infected_time = 0
    is_asymptomatic = 0 
    
    recover_time = 0 
    source = 'Imported' 
    class_type = '' 
    
    r = random.random()
    
    if r <= 0.1 :
        class_type = '1'
        isolation_time = today
        sym_time = today
        infected_time = today
        
    else:
        class_type = '5'
        isolation_time  = confirmed_time - confirmed_delay()
        sym_time = isolation_time - distr_logistic()
        infected_time = sym_time - import_infect2sym()
    
    recover_time = confirmed_time+ recovered_time() 
        
        
    case = [index, parent_index, infected_time, sym_time, isolation_time, confirmed_time, recover_time, source, class_type,is_asymptomatic]
    casesList.append(case)
    create_new_case(case)
    
    

def main(importList = []):

    for day in range(len(importList)):
        if day > last_day:
            return
        for import_num in range(importList[day]):
            imported_case(day)
           
    

if __name__ == "__main__":

    resultList = np.empty((0,last_day),int)
    dailyList = np.empty((0,last_day),int)
    for l in range(sim_num):
        if l%(sim_num/10) == 0:
            print("running...  "+ str(l)+' have done.')
        casesList = []
        
        index = 0
        
        main([1,2,1,0,1,2,3,3,3,2,
              0,0,2,1,0,0,0,0,1,0,
              0,0,0,0,0,1,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,
              0,1,1,2,2,1,3,1,9,5,
              9,9,9,10,18,35,25,29,39,19,
              48,32,40,29,22,42,23,10,16,19,
              8,9,7,5,3,3,4,3,1,0,
              0,0,0,0,0,1,1,0,0,0,
              0,0,0,0,2,0,0,0,0,0,
              0,0,0,0,0,0,0,0,1,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,1,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,1,4,1,3])

        name=['index', 'parent_index', 'infected_time', 'sym_time', 'isolation_time', 'confirmed_time', 'recover_time', 'source', 'class_type','is_asymptomatic']
        result = pd.DataFrame(columns=name,data=casesList)
        
        pre_y = [0]*last_day
        daily_pre_y = [0]*last_day
        for i in range(last_day):
            pre_y[i] = len(result[(result['confirmed_time']!=-100) & (result['confirmed_time']<=i) ])#& (result['source']=='Local')  & (result['recover_time']>i)
            daily_pre_y[i] = len(result[(result['confirmed_time']!=-100) & (result['confirmed_time']==i) ])#& (result['source']=='Local')  & (result['recover_time']>i)

        resultList = np.append(resultList, np.array([pre_y]), axis=0)
        dailyList = np.append(dailyList, np.array([daily_pre_y]), axis=0)

        '''
        # The effective reproduction in stages 2 and 3
        stage2 = result[(result['infected_time']>54) & (result['infected_time']<=76)]

        stage2_er1 = len(stage2[stage2.apply(lambda x: any(stage2['index'] == x.parent_index), axis=1)])
        stage2_er2 = len(stage2.groupby('parent_index'))
        stage2_er = stage2_er1/stage2_er2
        
        stage3 = result[(result['infected_time']>100) & (result['infected_time']<=132)]
        stage3_er1 = len(stage3[stage3.apply(lambda x: any(stage3['index'] == x.parent_index), axis=1)])
        stage3_er2 = len(stage3.groupby('parent_index'))
        stage3_er = stage3_er1/stage3_er2
        
        print(stage2_er,stage3_er)
        '''
    #Ground truth of the number of daily new confirmed cases   
    real_y = [1,3,4,4,5,7,10,13,16,18,18,18,24,28,30,33,
              40,42,44,45,48,56,65,70,73,75,79,82,83,84,
              87,87,88,89,91,94,96,100,104,106,108,110,115,128,136,148,
              158,164,176,185,198,210,224,240,264,311,343,383,430,453,
              507,556,629,680,728,797,837,867,908,969,1010,1053,
              1101,1161,1193,1241,1309,1367,1422,1460,1502,1549,
              1602,1641,1681,1720,1752,1789,1818,1859,1880,1920,
              1946,1960,1986,2007,2020,2038,2054,2064,2071,2086,
              2091,2101,2117,2127,2140,2150,2160,2163,2167,2171,
              2173,2175,2183,2192,2194,2195,2203,2216,2220,2231,
              2235,2241,2243,2247,2249,2258,2263,2265,2265,2269,
              2276,2291,2302,2309,2323,2325,2331,2338,2344,2362,
              2367,2377,2380,2382,2387,2391,2392,2394,2404,2405,
              2409,2416,2421,2427,2438,2450,2457,2463,
              2477,2488,2502,2512,2536,2563,2585,2597,2619,2635,2660,2662,2678,2687,
              2707,2722,2734,2749,2762,2775,2787,2801,2814,2816,2827,2836,2857,2864,
              2872,2878,2884
              ]

    
    real_y = real_y[0:min(len(real_y),last_day)]
    pre_cmu_mean = resultList.mean(0)[0:len(real_y)]
    
    real_y_daily_cases = [0]
    pre_y_daily_cases = [0]
    
    for i in range(1,len(real_y)):
        real_y_daily_cases.append(real_y[i]-real_y[i-1])
        pre_y_daily_cases.append(pre_cmu_mean[i]-pre_cmu_mean[i-1])
    
    error_daily_cases = np.mean(abs(np.array(pre_y_daily_cases)-np.array(real_y_daily_cases)))
    
    print('Daily new confirmed cases:')
    plt.plot(list(range(0,last_day)),dailyList.mean(0))
    plt.plot(list(range(0,len(real_y))),real_y_daily_cases)
    plt.show()
    
    print('Cumulative number of confirmed cases:')
    plt.plot(list(range(0,last_day)),resultList.mean(0))
    plt.plot(list(range(0,len(real_y))),real_y)
    plt.show()
    
    
    error = abs(pre_cmu_mean - real_y)
    print('Cumulative number error:')
    plt.plot(list(range(len(error))),error)
    plt.show()
    
    error = np.mean(error)
    
    print('daily error:'+ str(error_daily_cases))
    print('Cumulative error:'+str(error))
    
    error = abs(resultList[:,0:len(real_y)] - real_y).mean(0)
    error = np.mean(error)
