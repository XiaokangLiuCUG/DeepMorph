# Created by Xiaokang Liu on 24/07/2024 => xkliu@cug.edu.cn

import numpy as np
import matplotlib.pyplot as plt
    
def plot_results(coordinates,time, pattern,growth_rate,num_species):
    plt.figure(figsize=(7,7))  
    plt.clf()  
    ax = plt.gca() 
    ax.scatter(coordinates[:,0], coordinates[:,1],s = 81)
    ax.scatter(np.mean(coordinates[:,0]), np.mean(coordinates[:,1]), c='red',s = 100)
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    ax.set_aspect('equal') 
    plt.title("p = %s, time = %d,growth_rate= %.1f,num_sp = %d" % (pattern,time,growth_rate,num_species))
    plt.rcParams['font.family'] = 'Arial'  
    plt.savefig("morphological occupation, pattern = %s, time = %d,growth_rate= %.1f,num_sp = %d.pdf" % (pattern,time,growth_rate,num_species),format="svg")
    plt.show()  

def distances_sort(coordinates,point):
    distances = ((coordinates[:,0]-point[0])**2+(coordinates[:,1]-point[1])**2)**0.5
    sorted_index = sorted(range(len(distances)), key=lambda k: distances[k])
    sorted_coordinates = coordinates[sorted_index]
    return sorted_index 

def split_ex_rate(total_rate,num_taxa,selectivity):
    split_rate = np.array([])
#     for i in range(num_taxa):
#         split_rate = np.append(split_rate,2*(i+1)*total_rate /(num_taxa*(num_taxa+1)))
    split_rate = np.linspace(1/(100*(1+selectivity)), selectivity/(100*(1+selectivity)), 200,endpoint=True)
    return split_rate

def survivor_taxa(coordinates,sorted_index, extinction_rate,selectivity = 1,num_species = 200):
    m = selectivity
    #s = 1 / ((m-1)*num_species*extinction_rate+num_species)
#     higher_extinction_rate0 = np.repeat(m*s/2,int(200*extinction_rate))
#     higher_extinction_rate1 = split_ex_rate(m*s/2*int(200*extinction_rate),int(200*extinction_rate))
#     higher_extinction_rate = higher_extinction_rate0 + higher_extinction_rate1
    
#     higher_extinction_rate = np.repeat(m*s,int(200*extinction_rate))
    
#     lower_extinction_rate = np.repeat(s,round(200*(1-extinction_rate)))
#     extinction_rate_list = np.hstack((lower_extinction_rate,higher_extinction_rate))
    extinction_rate_list = split_ex_rate(1,200,selectivity)
    extinction_taxa_index = np.random.choice(sorted_index,int(200*extinction_rate),replace = False,p=extinction_rate_list.ravel())
    survivor_taxa_index = np.setdiff1d(sorted_index, extinction_taxa_index)
    survivor_taxa_coo = coordinates[survivor_taxa_index]
    return survivor_taxa_coo

def fit_multivariate_normal(coordinates):
    mean_est = np.mean(coordinates, axis=0)
    cov_est = np.cov(coordinates, rowvar=False)
    return mean_est, cov_est

def disparity(survivor_taxa_coo):
    sum_of_ranges = (np.max(survivor_taxa_coo[:,0]) - np.min(survivor_taxa_coo[:,0])) + (np.max(survivor_taxa_coo[:,1]) - np.min(survivor_taxa_coo[:,1]))
    sum_of_variances = np.var(survivor_taxa_coo[:,0]) + np.var(survivor_taxa_coo[:,1])
    sum_of_centroid_dis = np.sum(((survivor_taxa_coo[:,0] - np.mean(survivor_taxa_coo[:,0]))**2 + (survivor_taxa_coo[:,1] - np.mean(survivor_taxa_coo[:,1]))**2)**0.5) / survivor_taxa_coo.shape[0] 
    return [sum_of_ranges, sum_of_variances, sum_of_centroid_dis]
    
# def disparity_plot(total_disparity,extintion_rate):
#     marginal_disparity = total_disparity[:,0]
#     lateral_disparity = total_disparity[:,1]
#     random_survivor = total_disparity[:,2]
    
#     fig, ax = plt.subplots(1, 1,figsize=(7,7))
#     plt.plot(extintion_rate,marginal_disparity[:,1])
#     #plt.scatter(landmarks[:,0],landmarks[:,1],c = 'red')
#     plt.fill_between(extintion_rate,marginal_disparity[:,0],marginal_disparity[:,2],alpha=0.3,color='blue')
#     plt.xlim(0,1.0)
#     ax.invert_yaxis()
#     #plt.axis("equal")
#     plt.show()

def disparity_plot_one(total_disparity,extintion_rate):
    range_disparity = total_disparity[:,:,1]
    variance_disparity = total_disparity[:,:,1]
    
    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    plt.plot(extintion_rate,1-range_disparity[:,1])
    #plt.scatter(landmarks[:,0],landmarks[:,1],c = 'red')
    plt.fill_between(extintion_rate,1-range_disparity[:,0],1-range_disparity[:,2],alpha=0.3,color='blue')
    plt.xlim(0,1.0)
    #ax.invert_yaxis()
    #plt.axis("equal")
    plt.show()

def disparity_plot_four(non_se,two_se,five_se,ten_se,extintion_rate,text_label):

    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    plt.plot(extintion_rate,1-non_se[:,1],color='grey')
    plt.fill_between(extintion_rate,1-non_se[:,0],1-non_se[:,2],alpha=0.3,color='grey')
    
    plt.plot(extintion_rate,1-two_se[:,1],color=color[0])
    plt.fill_between(extintion_rate,1-two_se[:,0],1-two_se[:,2],alpha=0.3,color=color[0])
    
    plt.plot(extintion_rate,1-five_se[:,1],color=color[1])
    plt.fill_between(extintion_rate,1-five_se[:,0],1-five_se[:,2],alpha=0.3,color=color[1])
    
    plt.plot(extintion_rate,1-ten_se[:,1],color=color[2])
    plt.fill_between(extintion_rate,1-ten_se[:,0],1-ten_se[:,2],alpha=0.3,color=color[2])
    plt.xlim(0,1.0)
    #ax.invert_yaxis()
    #plt.axis("equal")
    plt.savefig("x-diversity loss,y-disparity,loss, type = %s.pdf" % (text_label),format="svg")
    plt.show()

def disparity_plot_all(total_disparity,extintion_rate,selectivity,c):
#     marginal_disparity = total_disparity[:,0]
#     lateral_disparity = total_disparity[:,1]
#     random_survivor = total_disparity[:,2]
    fig, axs = plt.subplots(3, 3,figsize=(12,12))
    fig.suptitle('simulations of selectivity = %s blue selective grey random'% str(selectivity-1))
    column = ['marginal_extinction','lateral_extinction','random_extinction']
    row = ['Soranges', 'Sovariances', 'centroid_dis']
    for i in range(3):
        for j in range(3):
                axs[i][j].plot(extintion_rate,1-total_disparity[i][:,j,1],color=c)
                axs[i][j].fill_between(extintion_rate,1-total_disparity[i][:,j,0],1-total_disparity[i][:,j,2],alpha=0.3,color=c)
                
                axs[i][j].plot(extintion_rate,1-total_disparity[i][:,2,1],color='grey')
                axs[i][j].fill_between(extintion_rate,1-total_disparity[i][:,2,0],1-total_disparity[i][:,2,2],alpha=0.3,color='grey')
                
                axs[i][j].set_title("%s of %s"%(column[j],row[i]))
#                 axs[i][j].set_xlim(0,1.0)
#                 axs[i][j].set_ylim(-0.4,1.0)
                #axs[i][j].invert_xaxis()
#     plt.plot(extintion_rate,marginal_disparity[:,1])
#     #plt.scatter(landmarks[:,0],landmarks[:,1],c = 'red')
#     plt.fill_between(extintion_rate,marginal_disparity[:,0],marginal_disparity[:,2],alpha=0.3,color='blue')
#     plt.axis("equal")
#     axs.set_xlim(0,1.0)
#     axs.set_ylim(-0.4,1.0)
    #axs.invert_yaxis()
    plt.savefig("x-diversity loss,y-disparity,loss, selectivity = %s.pdf" % (selectivity-1),format="svg")
    plt.show()

def selective_disparity(selectivity):
    total_disparity = np.array([])
    for n in extinction_rate_list:
        marginal_disparity = np.array([])
        lateral_disparity = np.array([])
        random_disparity = np.array([])
        for i in range(num_simulations):
            marginal_survivor_taxa_coo = survivor_taxa(np.asarray(marginal_sorted_index), extinction_rate = n, selectivity = selectivity,num_species=num_species)
            #lateral_survivor_taxa_coo = survivor_taxa(np.asarray(lateral_sorted_index), extinction_rate = n, selectivity = selectivity,num_species=num_species)
            marginal_disparity = np.append(marginal_disparity, disparity(marginal_survivor_taxa_coo))
            #lateral_disparity = np.append(lateral_disparity, disparity(lateral_survivor_taxa_coo))
        marginal_disparity = marginal_disparity.reshape(num_simulations,3)#[ranges, variances, centroid_dis]
        
        quantile_disparity = np.percentile(marginal_disparity, [2.5, 50, 97.5],axis=0)
        total_disparity = np.append(total_disparity,quantile_disparity)
    total_disparity = total_disparity.reshape(len(extinction_rate_list),3,3)
    return total_disparity / pre_extinction_disparity

def disparity_plot(three_folds,diversity,disparity_type,recovery_type,itime,growth_rate):
    time = list(range(11))
    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    plt.plot(time,three_folds[:,1],color='blue')
    plt.fill_between(time,three_folds[:,0],three_folds[:,2],alpha=0.3,color='blue')
    
    plt.plot(time,diversity[:,1],color='red')
    plt.fill_between(time,diversity[:,0],diversity[:,2],alpha=0.3,color='grey')
    
    plt.title("disparity_type = %s, recovery_type = %s, time = %d" % (disparity_type,recovery_type,itime))
    plt.savefig("x-time,y-disparity,increase, disparity_type = %s, recovery_type = %s, time = %d, growth_rate = %f.pdf" % (disparity_type,recovery_type,itime,growth_rate),format="svg")
    plt.show()
    
np.random.seed(3)#3
num_species = 200
start_Stdeviation = 0.5
coordinates0 = np.random.normal(0, start_Stdeviation,[num_species,2])# dataset generation
selectivity = 100 # 0.9 # 2/3
lateral_index = 1.5
net_diversification_rate = 0.1
num_simulations = 1000
print(coordinates0.shape)
extinction_rate = 0.8 #[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]#
N0 = int(num_species * (1-extinction_rate)) # 40
plot_results(coordinates0,0,'pre-extinction',growth_rate=1,num_species=200)

ext_rate0 = extinction_rate + np.random.uniform(10, 15)/200
marginal_sorted_index0 = distances_sort(coordinates0,np.array([np.mean(coordinates0[:,0]),np.mean(coordinates0[:,1])]))
marginal_survivor_taxa_coo0 = survivor_taxa(coordinates0,np.asarray(marginal_sorted_index0), extinction_rate = ext_rate0, selectivity = selectivity,num_species=num_species)
lateral_sorted_index0 = distances_sort(coordinates0,np.array([np.min(coordinates0[:,0])*lateral_index,np.mean(coordinates0[:,1])]))
lateral_survivor_taxa_coo0 = survivor_taxa(coordinates0,np.asarray(lateral_sorted_index0), extinction_rate = ext_rate0, selectivity = selectivity,num_species=num_species)
plot_results(marginal_survivor_taxa_coo0,0,'pre-ext marg_survivor_taxa_coo0',growth_rate=1,num_species = int(ext_rate0*200))
plot_results(lateral_survivor_taxa_coo0,0,'pre-ext late_survivor_taxa_coo0',growth_rate=1,num_species = int(ext_rate0*200))
print("lateral_survivor_taxa_coo0 average centroid：",np.mean(lateral_survivor_taxa_coo0,axis=0))
pre_extinction_disparity = np.array([])
pre_ext_diversity = np.array([])
lateral_disparity = marginal_disparity = np.array([])
pre_marginal_ext_disparity = pre_lateral_ext_disparity = np.array([])
#total_disparity = np.append(total_disparity,pre_extinction_disparity)
for i in range(num_simulations):
    coordinates = np.random.normal(0, start_Stdeviation,[num_species,2])
    pre_extinction_disparity = np.append(pre_extinction_disparity, disparity(coordinates))
    ext_rate = extinction_rate + np.random.normal(0, 0.05)
    pre_ext_diversity = np.append(pre_ext_diversity,200-int(ext_rate*200))
    
    marginal_sorted_index = distances_sort(coordinates,np.array([np.mean(coordinates[:,0]),np.mean(coordinates[:,1])]))
    marginal_survivor_taxa_coo = survivor_taxa(coordinates,np.asarray(marginal_sorted_index), extinction_rate = ext_rate, selectivity = selectivity,num_species=num_species)
    marginal_disparity = np.append(marginal_disparity, disparity(marginal_survivor_taxa_coo))
    
    lateral_sorted_index= distances_sort(coordinates,np.array([np.min(coordinates[:,0])*lateral_index,np.mean(coordinates[:,1])]))
    lateral_survivor_taxa_coo = survivor_taxa(coordinates,np.asarray(lateral_sorted_index), extinction_rate = ext_rate, selectivity = selectivity,num_species=num_species)
    lateral_disparity = np.append(lateral_disparity, disparity(lateral_survivor_taxa_coo))
    
#plot_results(marginal_survivor_taxa_coo,1,"marginal_survivor_taxa_coo")
#plot_results(lateral_survivor_taxa_coo,1,"lateral_survivor_taxa_coo")

pre_extinction_disparity = pre_extinction_disparity.reshape(num_simulations,3)
marginal_disparity = marginal_disparity.reshape(num_simulations,3)
lateral_disparity = lateral_disparity.reshape(num_simulations,3)

quantile_marginal_disparity = np.percentile(marginal_disparity, [2.5, 50, 97.5],axis=0)
quantile_lateral_disparity = np.percentile(lateral_disparity, [2.5, 50, 97.5],axis=0)
quantile_pre_extinction_disparity = np.percentile(pre_extinction_disparity, [2.5, 50, 97.5],axis=0)
quantile_pre_ext_diversity = np.percentile(pre_ext_diversity, [2.5, 50, 97.5],axis=0)
start_diversity = np.percentile(200 * (1+np.random.normal(0, 0.05,3000)), [2.5, 50, 97.5])

pre_marginal_ext_disparity = np.append(quantile_pre_extinction_disparity,quantile_marginal_disparity)
pre_lateral_ext_disparity = np.append(quantile_pre_extinction_disparity,quantile_lateral_disparity)
pre_ext_diversity = np.append(start_diversity,quantile_pre_ext_diversity)

expansion_factors = [1,2,0.5]#backfill, expansion, and stagnation mode
metrics = ["range","variance","centroid_distance"]
recovery_patterns = ["refilling mode", "expansion mode","stagnation mode"]


def recovery_pattern(diversity_ratio =1, growth_rate =1, Tn_centroid = [0,0], shift_ratio = 0, expansion_factor = 1,pre_ext_disparity=pre_marginal_ext_disparity,survivor_taxa_coo=marginal_survivor_taxa_coo,recovery_pattern="backfill mode"):
    total_disparity = pre_ext_disparity
    total_diversity = start_diversity
    for i in range(0,10):
        T = i
        Tn_centroid = shift_ratio * (1 / 9)* (9-T)#[shift_ratio * (1.2 / 9)* (9-T),0]
        
        Tn_Stdeviation = np.mean(np.std(survivor_taxa_coo, axis=0)) + expansion_factor * (start_Stdeviation-np.mean(np.std(survivor_taxa_coo, axis=0))) / (1 + ((num_species - N0) / N0) * np.exp(-growth_rate * T))
        sim_disparities = np.array([])
        sim_diversities = np.array([])
        for m in range(num_simulations):
            Tn_species = diversity_ratio * num_species / (1 + ((num_species - N0) / N0) * np.exp(-growth_rate * T)) + np.random.uniform(-8 , 8)#np.random.uniform(-T * 2 , T * 2)
            #Tn_species = diversity_ratio * (1+np.random.normal(0, 0.08)) * (num_species * net_diversification_rate * T) #add some randomness 
            T_species = np.random.normal(Tn_centroid, Tn_Stdeviation,[int(Tn_species),2])
            T_disparity = disparity(T_species)
            sim_disparities = np.append(sim_disparities, T_disparity)
            sim_diversities = np.append(sim_diversities, Tn_species)
        sim_disparities = sim_disparities.reshape(num_simulations,3)
        quantile_diversity = np.percentile(sim_diversities, [2.5, 50, 97.5],axis=0)
        quantile_disparity = np.percentile(sim_disparities, [2.5, 50, 97.5],axis=0)
        total_disparity = np.append(total_disparity,quantile_disparity)
        total_diversity = np.append(total_diversity,quantile_diversity)
        print("Number of species：Time：",Tn_species, T)
        if i in [0,4,9]:
            print("time:",i)
            print(T_species.shape)
            plot_results(T_species,i,recovery_pattern,growth_rate,Tn_species)
    
    total_disparity = total_disparity.reshape(11,3,3) / np.mean(pre_extinction_disparity,axis=0)
    total_diversity = total_diversity.reshape(11,3) / 200
    print(total_diversity[:,1])
    for i in range(2):
        disparity_plot(total_disparity[:,:,i],total_diversity,metrics[i],recovery_pattern,i+4,diversity_ratio)

# backfill mode, lateral_marginal_disparity
recovery_pattern(diversity_ratio = np.random.uniform(0.8, 1),growth_rate = 1, shift_ratio = np.mean(lateral_survivor_taxa_coo0,axis=0),expansion_factor = expansion_factors[0],survivor_taxa_coo=lateral_survivor_taxa_coo,pre_ext_disparity=quantile_pre_extinction_disparity,recovery_pattern=recovery_patterns[0])
# expansion mode
recovery_pattern(diversity_ratio = np.random.uniform(0.7, 0.8),growth_rate = 0.5,shift_ratio = 0,expansion_factor = expansion_factors[1],survivor_taxa_coo=marginal_survivor_taxa_coo,pre_ext_disparity=quantile_pre_extinction_disparity,recovery_pattern=recovery_patterns[1])
# stagnation mode
recovery_pattern(diversity_ratio = np.random.uniform(0.6, 0.7),growth_rate = 1.3,shift_ratio = 0,expansion_factor = expansion_factors[2],survivor_taxa_coo=marginal_survivor_taxa_coo,pre_ext_disparity=quantile_pre_extinction_disparity,recovery_pattern=recovery_patterns[2])
