import numpy as np
import matplotlib.pyplot as plt
#import random
import numpy as np
import matplotlib.pyplot as plt
#import random
def plot_results(coordinates,lateral_index):
    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    #ax.invert_yaxis()
    plt.scatter(coordinates[:,0],coordinates[:,1])
    plt.scatter(np.mean(coordinates[:,0]),np.mean(coordinates[:,1]),c = 'red')
    plt.scatter(np.max(coordinates[:,0])*lateral_index,np.mean(coordinates[:,1]),c = 'blue')
    plt.axis("equal")
    plt.show()

def distances_sort(coordinates,point):
    distances = ((coordinates[:,0]-point[0])**2+(coordinates[:,1]-point[1])**2)**0.5
    sorted_index = sorted(range(len(distances)), key=lambda k: distances[k])
    sorted_coordinates = coordinates[sorted_index]
    return sorted_index #距离从小到大排序

def split_ex_rate(total_rate,num_taxa,selectivity):
    split_rate = np.array([])
#     for i in range(num_taxa):
#         split_rate = np.append(split_rate,2*(i+1)*total_rate /(num_taxa*(num_taxa+1)))
    split_rate = np.linspace(1/(100*(1+selectivity)), selectivity/(100*(1+selectivity)), 200,endpoint=True)
    return split_rate

def survivor_taxa(sorted_index, extinction_rate,selectivity = 1,num_species = 200):
    m = selectivity# 倍数
    #s = 1 / ((m-1)*num_species*extinction_rate+num_species)
#     higher_extinction_rate0 = np.repeat(m*s/2,int(200*extinction_rate))
#     higher_extinction_rate1 = split_ex_rate(m*s/2*int(200*extinction_rate),int(200*extinction_rate))
#     higher_extinction_rate = higher_extinction_rate0 + higher_extinction_rate1
#     higher_extinction_rate = np.repeat(m*s,int(200*extinction_rate)) 
#     lower_extinction_rate = np.repeat(s,round(200*(1-extinction_rate)))
#     extinction_rate_list = np.hstack((lower_extinction_rate,higher_extinction_rate))
    extinction_rate_list = split_ex_rate(1,200,selectivity)
    extinction_taxa_index = np.random.choice(sorted_index,int(200*extinction_rate),False,p=extinction_rate_list.ravel())
    survivor_taxa_index = np.setdiff1d(sorted_index, extinction_taxa_index)
    survivor_taxa_coo = coordinates[survivor_taxa_index]
    return survivor_taxa_coo

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
    
np.random.seed(3)
num_species = 200
coordinates = np.random.normal(0, 0.5,[num_species,2])# dataset generation
selectivities = [2,6,11]# 0.9 # 2/3
lateral_index = 0.6
num_simulations = 1000
print(coordinates.shape)
extinction_rate_list = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]#
plot_results(coordinates,lateral_index)
pre_extinction_disparity = disparity(coordinates)#[ranges, variances, centroid_dis]
marginal_sorted_index = distances_sort(coordinates,np.array([np.mean(coordinates[:,0]),np.mean(coordinates[:,1])]))
lateral_sorted_index = distances_sort(coordinates,np.array([np.max(coordinates[:,0])*lateral_index,np.mean(coordinates[:,1])]))
color = ['orange',"blue","red"]

def disparity_plot_x(total_diversity_rate,disparity_loss_list,selectivity,c):
    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    #plt.plot(disparity_loss_list,total_diversity_rate[0][:,0],color='red')
    plt.plot(disparity_loss_list,total_diversity_rate[0][:,1],color=c)
    plt.fill_between(disparity_loss_list,total_diversity_rate[0][:,0],total_diversity_rate[0][:,2],alpha=0.3,color=c)
    
    plt.plot(disparity_loss_list,total_diversity_rate[1][:,1],color='grey')
    plt.fill_between(disparity_loss_list,total_diversity_rate[1][:,0],total_diversity_rate[1][:,2],alpha=0.3,color='grey')
    plt.title("Range vars x-dis loss,y-di loss selectivity = %s" % (selectivity-1))
    #axs.invert_yaxis()
    plt.axis("equal")
    ax.set_xlim(0,1.0)
    ax.set_ylim(0,1.0)
    plt.savefig("Range vars x-dis loss,y-div loss selectivity = %s.pdf" % (selectivity-1),format="svg")
    plt.show()

def disparity_plot_all(total_disparity,extintion_rate,extinction_rate):
#     marginal_disparity = total_disparity[:,0]
#     lateral_disparity = total_disparity[:,1]
#     random_survivor = total_disparity[:,2]
    fig, axs = plt.subplots(3, 3,figsize=(15,15))
    fig.suptitle('simulations of extinction rate = %s blue selective red random'% str(extinction_rate))
    column = ['marginal_extinction','lateral_extinction','random_extinction']
    row = ['SOranges', 'SOvariances', 'centroid_dis']
    for i in range(3):
        for j in range(3):
                axs[i][j].plot(selectivities,total_disparity[i][:,j,1]-1,color='blue')
                axs[i][j].fill_between(selectivities,total_disparity[i][:,j,0]-1,total_disparity[i][:,j,2]-1,alpha=0.3,color='blue')
                
                axs[i][j].plot(selectivities,total_disparity[i][:,2,1]-1,color='red')
                axs[i][j].fill_between(selectivities,total_disparity[i][:,2,0]-1,total_disparity[i][:,2,2]-1,alpha=0.3,color='red')
                
                axs[i][j].set_title("%s of %s"%(column[j],row[i]))
                #axs[i][j].invert_yaxis()
    plt.xlim(0,1.0)
    #axs.invert_yaxis()
    #plt.axis("equal")
    plt.show()

def div_vs_dis(selectivities,disparity_loss_list,color)
    for selectivity in selectivities:
        total_ranges = np.array([])
        total_variances = np.array([])
        total_centroid_dis = np.array([])
        for n in extinction_rate_list:
            marginal_disparity = np.array([])
            lateral_disparity = np.array([])
            random_disparity = np.array([])
            for i in range(num_simulations):
                marginal_survivor_taxa_coo = survivor_taxa(np.asarray(marginal_sorted_index), extinction_rate = n, selectivity = selectivity,num_species=num_species)
                lateral_survivor_taxa_coo = survivor_taxa(np.asarray(lateral_sorted_index), extinction_rate = n, selectivity = selectivity,num_species=num_species)
                marginal_disparity = np.append(marginal_disparity, disparity(marginal_survivor_taxa_coo))
                lateral_disparity = np.append(lateral_disparity, disparity(lateral_survivor_taxa_coo))
                random_survivor = np.random.choice(np.arange(0,200),round(200*(1-n)),False)
                random_disparity = np.append(random_disparity, disparity(coordinates[random_survivor]))
            marginal_disparity = marginal_disparity.reshape(num_simulations,3)#[ranges, variances, centroid_dis]
            lateral_disparity = lateral_disparity.reshape(num_simulations,3)
            random_disparity = random_disparity.reshape(num_simulations,3)
            list1 = [marginal_disparity,lateral_disparity,random_disparity]
            for m in range(3):
                total_ranges = np.append(total_ranges, np.percentile(list1[m][:,0], [2.5, 50, 97.5]))
                total_variances = np.append(total_variances, np.percentile(list1[m][:,1], [2.5, 50, 97.5]))
                total_centroid_dis = np.append(total_centroid_dis, np.percentile(list1[m][:,2], [2.5, 50, 97.5]))
        total_ranges = total_ranges.reshape(len(extinction_rate_list),3,3) / pre_extinction_disparity[0]
        total_variances = total_variances.reshape(len(extinction_rate_list),3,3)/ pre_extinction_disparity[1]
        total_centroid_dis = total_centroid_dis.reshape(len(extinction_rate_list),3,3) / pre_extinction_disparity[2]
        total_disparity = [total_ranges, total_variances, total_centroid_dis]
        #disparity_plot(total_variances,extinction_rate_list)# 三种灭绝模式，三种计算方式
        disparity_plot_all(total_disparity,extinction_rate_list,selectivity,color[selectivity])#total_ranges total_centroid_dis total_variances
        
def dis_vs_div(selectivities,disparity_loss_list)
    for selectivity in range(len(selectivities)):
        marginal_rate = np.array([])
        random_rate = np.array([])
        for x in disparity_loss_list:
            total_ranges = np.array([])
            total_variances = np.array([])
            total_centroid_dis = np.array([])
            var_diversity_loss = np.array([])
            random_var_diversity_loss = np.array([])
            distance_diversity_loss = np.array([])
            for n in range(2,199):
                marginal_disparity = np.array([])
                lateral_disparity = np.array([])
                random_disparity = np.array([])
                for i in range(num_simulations):
                    marginal_survivor_taxa_coo = survivor_taxa(np.asarray(marginal_sorted_index), extinction_rate = n/num_species, selectivity = selectivities[selectivity],num_species=num_species)
                    #lateral_survivor_taxa_coo = survivor_taxa(np.asarray(lateral_sorted_index), extinction_rate = n, selectivity = selectivity,num_species=num_species)
                    marginal_disparity = np.append(marginal_disparity, disparity(marginal_survivor_taxa_coo))
                    #lateral_disparity = np.append(lateral_disparity, disparity(lateral_survivor_taxa_coo))
                    random_survivor = np.random.choice(np.arange(0,200),round(200-n),False)
                    random_disparity = np.append(random_disparity, disparity(coordinates[random_survivor]))
                #print("marginal_disparity",np.percentile(marginal_disparity, [2.5])/pre_extinction_disparity )

                x1 = np.percentile(marginal_disparity, [2.5])/pre_extinction_disparity
                x2 = np.percentile(marginal_disparity, [97.5])/pre_extinction_disparity
                x3 = np.percentile(random_disparity, [2.5])/pre_extinction_disparity
                x4 = np.percentile(random_disparity, [97.5])/pre_extinction_disparity
                if (x1[0] < (1-x)) and ((1-x) < x2[0]):
                    var_diversity_loss = np.append(var_diversity_loss,n)
                if (x3[0] < (1-x)) and ((1-x) < x4[0]):
                    random_var_diversity_loss = np.append(random_var_diversity_loss,n)
            print("disparity loss rate %s "% str(x))
            marginal_rate = np.append(marginal_rate, np.percentile(var_diversity_loss/200, [2.5, 50, 97.5]))
            random_rate = np.append(random_rate, np.percentile(random_var_diversity_loss/200, [2.5, 50, 97.5]))

        # plot
        marginal_rate = marginal_rate.reshape(len(disparity_loss_list),3)
        random_rate = random_rate.reshape(len(disparity_loss_list),3)
        total_diversity_rate = [marginal_rate, random_rate]
        disparity_plot_x(total_diversity_rate,disparity_loss_list,selectivities[selectivity],color[selectivity])

if __name__ == '__main__':
    # parameters
    np.random.seed(3)
    num_species = 200
    coordinates = np.random.normal(0, 0.5,[num_species,2])# dataset generation
    selectivity = 0.6
    lateral_index = 0.6
    num_simulations = 500
    print(coordinates.shape)
    extinction_rate_list = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]#
    disparity_loss_list = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
    selectivities = [2,5,10]
    plot_results(coordinates,lateral_index)
    pre_extinction_disparity = disparity(coordinates)#[ranges, variances, centroid_dis]
    marginal_sorted_index = distances_sort(coordinates,np.array([np.mean(coordinates[:,0]),np.mean(coordinates[:,1])]))
    lateral_sorted_index = distances_sort(coordinates,np.array([np.max(coordinates[:,0])*lateral_index,np.mean(coordinates[:,1])]))
    color = ['orange',"blue","red"]
    
    # diversity loss vs.  disparity loss
    div_vs_dis(selectivities,disparity_loss_list,color)
    
    # disparity loss vs. diversity loss
    dis_vs_div(selectivities,disparity_loss_list,color)