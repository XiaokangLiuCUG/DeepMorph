import os
import random
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt,atan2
import xlwt

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 255)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

# contours
def get_coordinates(path):
    image = Image.open(path).convert("L")
    img_array = np.array(image)
    img_array[img_array>=150]=255
    img_array[img_array<150]=0
    img_array = np.pad(img_array, 2, pad_with)
    icount = np.array([])
    jcount = np.array([])
    for i in range(img_array.shape[0]-1):
        #print(i,img_array[i])
        for j in range(img_array.shape[1]-1):
            if img_array[i][j] != img_array[i][j+1]:
                if img_array[i][j]==0:
                    icount = np.append(icount,i) 
                    jcount = np.append(jcount,j)
                else:
                    icount = np.append(icount,i) 
                    jcount = np.append(jcount,j+1)
            if img_array[i][j] != img_array[i+1][j]:
                if img_array[i][j]==0:
                    icount = np.append(icount,i) 
                    jcount = np.append(jcount,j)
                else:
                    icount = np.append(icount,i+1)
                    jcount = np.append(jcount,j)
    total_coordinates = np.stack((jcount,icount),axis = 1)
    coordinates = []# np.np.unique
    for i in total_coordinates.tolist():
        if i not in coordinates:
            coordinates.append(i)
    return np.array(coordinates,dtype=int)

###### four fixed landmarks ######
def four_landmarks(coordinates):
    X,Y = coordinates[:,0],coordinates[:,1]
    num_top_index = [i for i,x in enumerate(Y.tolist()) if x == int(np.min(Y))]
    middle_top_index = num_top_index[int(len(num_top_index) / 2)]
    top_point_x, top_point_y = X[middle_top_index],Y[middle_top_index]
    
    num_down_index = [i for i,x in enumerate(Y.tolist()) if x == int(np.max(Y))]
    middle_down_index = num_down_index[int(len(num_down_index) / 2)]
    down_point_x, down_point_y = X[middle_down_index],Y[middle_down_index]
    
    if not (abs(top_point_x-(np.min(X) + np.max(X))/2)) / ((np.min(X) + np.max(X))/2-np.min(X)) < 0.2:
        top_point_x = int((np.min(X) + np.max(X))/2)
        top_middel_index = [i for i,x in enumerate(X.tolist()) if x == top_point_x]
        if len(top_middel_index)>2:
            print("********* need to revise *********")
        top_point_y = Y[top_middel_index[0]]
        #down_point_x, down_point_y = top_point_x, Y[top_middel_index[-1]]
        
    if not (abs(down_point_x-(np.min(X) + np.max(X))/2))/ ((np.min(X) + np.max(X))/2-np.min(X)) < 0.2:
        down_point_x = int((np.min(X) + np.max(X))/2)
        down_middel_index = [i for i,x in enumerate(X.tolist()) if x == down_point_x]
        if len(down_middel_index)>2:
            print("********* need to revise *********")
        down_point_y = Y[down_middel_index[-1]]
        #down_point_x, down_point_y = top_point_x, Y[down_point_y]
        
    left_point_y = right_point_y = int((np.min(Y) + np.max(Y))/2) #int((top_point_y + down_point_y) / 2)#whole shell mean Y
    left_rigeht_point_index = [i for i,x in enumerate(Y.tolist()) if x == left_point_y]
    left_point_index = [i for i in left_rigeht_point_index if X[i] < top_point_x]
    right_point_index = [i for i in left_rigeht_point_index if X[i] > top_point_x]
    left_point_x = X[left_point_index[int((len(left_point_index)-1)/2)]]
    right_point_x = X[right_point_index[int((len(right_point_index)-1)/2)]]
    landmarks_x = [top_point_x,down_point_x,left_point_x,right_point_x]
    landmarks_y = [top_point_y,down_point_y,left_point_y,right_point_y]
    four_landmarks_points = np.stack((landmarks_x,landmarks_y),axis = 1)

    return four_landmarks_points


###### find left right landmarks ######
def left_right_landmarks(X,Y):
    num_left_index = [i for i,x in enumerate(X.tolist()) if x == int(np.min(X))]
    middle_left_index = num_left_index[int(len(num_left_index) / 2)]
    num_right_index = [i for i,x in enumerate(X.tolist()) if x == int(np.max(X))]
    middle_right_index = num_right_index[int(len(num_right_index) / 2)]
    #print("top_index",top_index)
    left_point_x, left_point_y = X[middle_left_index],Y[middle_left_index]
    right_point_x, right_point_y = X[middle_right_index],Y[middle_right_index]
    #print("左顶点：",left_point_x, left_point_y)

    landmarks_x = [left_point_x,right_point_x]
    landmarks_y = [left_point_y,right_point_y]
    two_landmarks_points = np.stack((landmarks_x,landmarks_y),axis = 1)
    return two_landmarks_points

###### find top down landmarks ######
def top_down_landmarks(coordinates):
    X,Y = coordinates[:,0],coordinates[:,1]
    num_top_index = [i for i,x in enumerate(Y.tolist()) if x == int(np.min(Y))]
    middle_top_index = num_top_index[int(len(num_top_index) / 2)]
    top_point_x, top_point_y = X[middle_top_index],Y[middle_top_index]
    
    num_down_index = [i for i,x in enumerate(Y.tolist()) if x == int(np.max(Y))]
    middle_down_index = num_down_index[int(len(num_down_index) / 2)]
    down_point_x, down_point_y = X[middle_down_index],Y[middle_down_index]
    
    if not (abs(top_point_x-(np.min(X) + np.max(X)))/2)/ ((np.min(X) + np.max(X))/2-np.min(X)) < 0.2:
        print("out of upper boundary")
        top_point_x = int((np.min(X) + np.max(X)/2))
        top_middel_index = [i for i,x in enumerate(X.tolist()) if x == top_point_x]
        if len(top_middel_index)>2:
            print("********* needs to revise *********")
        top_point_y = Y[top_middel_index[0]]
        #down_point_x, down_point_y = top_point_x, Y[top_middel_index[-1]]
    
    if not (abs(down_point_x-(np.min(X) + np.max(X)))/2)/ ((np.min(X) + np.max(X))/2-np.min(X)) < 0.2:
        down_point_x = int((np.min(X) + np.max(X)/2))
        down_middel_index = [i for i,x in enumerate(X.tolist()) if x == down_point_x]
        if len(down_middel_index)>2:
            print("********* needs to revise *********")
        down_point_y = Y[down_middel_index[-1]]
        down_point_x, down_point_y = top_point_x, Y[down_point_y[-1]]
        print("out of lower boundary")
    
    #print(top_point_x,top_point_y)
    landmarks_x = [top_point_x,down_point_x]
    landmarks_y = [top_point_y,down_point_y]
    two_landmarks_points = np.stack((landmarks_x,landmarks_y),axis = 1)
    return two_landmarks_points

def semi_landmarks(coordinates,num_semilandmarks):
    choose = np.array([],dtype=int)
    distance = coordinates.shape[0] / (num_semilandmarks+1)
    for i in range(num_semilandmarks):
        choose = np.append(choose,coordinates[round(distance * (1+i))],axis=0)
    choose = choose.reshape(int(choose.shape[0] / 2),2)
    return choose

def length(v):           # length calculation
    return sqrt(v[0]**2 + v[1]**2) 

def to_polar(vector, centered_x, centered_y):   #Cartesian coordinates to polar coordinates
    x, y = vector[0]-centered_x, vector[1]-centered_y
    angle = atan2(y,x)
    #return (length(vector), angle)
    return angle

def sort_coordinates(coordinates,centered_x, centered_y):
    angles = []
    for i in range(coordinates.shape[0]):
        angles.append(to_polar(coordinates[i], centered_x, centered_y))
    sorted_angles = sorted(range(len(angles)), key=lambda k: angles[k])
    sorted_coordinates = coordinates[sorted_angles]
    return sorted_coordinates

def convcave_judgement(coordinates):
    convcave = False
    for i in range(1,coordinates.shape[0]):
        if coordinates[i][0]-coordinates[i-1][0] >= 3:
            print("concave polygon")
            convcave=True
            break
    return convcave

def distance_sort(start):
    contain = [start[0].tolist()]
    for i in range(start.shape[0]):
        if i == 0:
            index = 0
            point = contain[i]
        else:
            point = start[index]
        distance = (start[:,0]-point[0])**2+(start[:,1]-point[1])**2
        for m in range(1,20):
            index = np.argsort(distance)[m]
            if start[index].tolist() not in contain:
                contain.append(start[index].tolist())
                point = start[index].tolist()
                break
    if len(contain) < start.shape[0]*0.9:
        print("************ polygon needs revise ************")
    return np.array(contain)

# dorsal_view four quadrant coordinates split
def coordinates_split(coordinates,num_semilandmarks,zone, num):
    quadrant = np.array([])
    X,Y = coordinates[:,0],coordinates[:,1]
    for i in range(len(X)):
        if num == 1:
            if (np.min(zone[:,0]) <= X[i]) and (Y[i] <= np.max(zone[:,1])):
                quadrant = np.append(quadrant,coordinates[i])
        if num == 2:
            if (X[i] <= np.max(zone[:,0])) and (Y[i] <= np.max(zone[:,1])):
                quadrant = np.append(quadrant,coordinates[i])
        if num == 3:
            if (X[i] <= np.max(zone[:,0])) and (np.min(zone[:,1]) <= Y[i]):
                quadrant = np.append(quadrant,coordinates[i])
        if num == 4:
            if (np.min(zone[:,0]) <= X[i]) and (np.min(zone[:,1]) <= Y[i]):
                quadrant = np.append(quadrant,coordinates[i])
    quadrant = quadrant.reshape(int(quadrant.shape[0] / 2),2)
    #print(len(quadrant_x))
    return semi_landmarks(quadrant,num_semilandmarks)

def dorsal_view(four_landmarks_points,coordinates,num_upper_semilandmarks,num_lower_semilandmarks):
    centered_x = four_landmarks_points[0][0]
    centered_y = four_landmarks_points[2][1]
    sorted_coordinates = sort_coordinates(coordinates,centered_x, centered_y)
    first_zone = np.array([four_landmarks_points[0],four_landmarks_points[3]])
    second_zone = np.array([four_landmarks_points[0],four_landmarks_points[2]])
    thrid_zone = np.array([four_landmarks_points[1],four_landmarks_points[2]])
    fourth_zone = np.array([four_landmarks_points[1],four_landmarks_points[3]])

    first_quadrant = coordinates_split(sorted_coordinates,num_upper_semilandmarks,first_zone,1)
    second_quadrant = coordinates_split(sorted_coordinates,num_upper_semilandmarks,second_zone,2)
    thrid_quadrant = coordinates_split(sorted_coordinates,num_lower_semilandmarks,thrid_zone,3)
    fourth_quadrant = coordinates_split(sorted_coordinates,num_lower_semilandmarks,fourth_zone,4)
    total_landmarks = np.concatenate(([four_landmarks_points[0]],first_quadrant,[four_landmarks_points[3]],fourth_quadrant,[four_landmarks_points[1]],thrid_quadrant,[four_landmarks_points[2]],second_quadrant))
    return total_landmarks

###### b view semi-landmarks ######
def left_right_coordinates_split(coordinates,num_semilandmarks,two_landmarks):
    index2 = coordinates.tolist().index(two_landmarks[1].tolist())
    quadrant1 = quadrant2 = np.array([],dtype=int)
    for i in range(coordinates.shape[0]):
        if i < index2:
            quadrant1 = np.append(quadrant1,coordinates[i])
        else:
            quadrant2 = np.append(quadrant2,coordinates[i])
    quadrant1 = quadrant1.reshape(int(quadrant1.shape[0] / 2),2)
    quadrant2 = quadrant2.reshape(int(quadrant2.shape[0] / 2),2)
    right_landmarks = semi_landmarks(quadrant1,num_semilandmarks)
    left_landmarks = semi_landmarks(quadrant2,num_semilandmarks)    
    return right_landmarks,left_landmarks

def front_view(two_landmarks_points,coordinates,num_front_semilandmarks):
    centered_x = int((two_landmarks_points[0][0]+two_landmarks_points[1][0])/2)
    centered_y = two_landmarks_points[0][1]
    #print(centered_x,centered_y)
    sorted_coordinates = sort_coordinates(coordinates,centered_x, centered_y)
    right_landmarks,left_landmarks = left_right_coordinates_split(sorted_coordinates,num_front_semilandmarks,two_landmarks_points)
    total_landmarks = np.concatenate(([two_landmarks_points[0]],right_landmarks,[two_landmarks_points[1]],left_landmarks))
    return total_landmarks

###### c view semi-landmarks ######
def top_down_coordinates_split(coordinates,num_semilandmarks,two_landmarks_points):
    index1 = coordinates.tolist().index(two_landmarks_points[0].tolist())
    index2 = coordinates.tolist().index(two_landmarks_points[1].tolist())
    quadrant1 = quadrant2 = np.array([],dtype=int)
    for i in range(index2,coordinates.shape[0]+index2):
        if i >= coordinates.shape[0]:
            i = i-coordinates.shape[0]
        if index2 > i > index1:
            quadrant2 = np.append(quadrant2,coordinates[i])
        else:
            quadrant1 = np.append(quadrant1,coordinates[i])

    quadrant1 = quadrant1.reshape(int(quadrant1.shape[0] / 2),2)
    quadrant2 = quadrant2.reshape(int(quadrant2.shape[0] / 2),2)
    left_landmarks = semi_landmarks(quadrant1,num_semilandmarks)
    right_landmarks = semi_landmarks(quadrant2,num_semilandmarks)    
    return right_landmarks,left_landmarks

def lateral_view(two_landmarks_points, coordinates,num_lateral_semilandmarks):
    centered_x = two_landmarks_points[0][0]
    centered_y = int((two_landmarks_points[1][1]+two_landmarks_points[0][1])/2)
    sorted_coordinates = sort_coordinates(coordinates,centered_x, centered_y)
    convcave = convcave_judgement(sorted_coordinates)
    if convcave:
        sorted_coordinates = distance_sort(sorted_coordinates)
        if to_polar(sorted_coordinates[50], centered_x, centered_y) > 0:
            sorted_coordinates = sorted_coordinates[::-1]
    right_landmarks,left_landmarks = top_down_coordinates_split(sorted_coordinates,num_lateral_semilandmarks,two_landmarks_points)
    for i in range(1,sorted_coordinates.shape[0]):
        if sorted_coordinates[i][0]-sorted_coordinates[i-1][0] >= 5:
            sorted_coordinates = distance_sort(coordinates)
            print("concave polygon")
            break
            right_landmarks,left_landmarks = top_down_coordinates_split_con(sorted_coordinates,num_lateral_semilandmarks,two_landmarks_points)
    total_landmarks = np.concatenate(([two_landmarks_points[0]],right_landmarks,[two_landmarks_points[1]],left_landmarks))
    return total_landmarks

def readpaths(path):
    with open(path, 'r', encoding='utf-8') as f:
        mylist = []
        for line in f.readlines():
            line = line.strip('\n')  
            mylist.append(line)
    return mylist

def get_list(path):
    img_lists = []
    species_names = []
    total_views = []
    img_names = []
    genus = []
    for root,dirs,files in os.walk(path):
        for file in files:
            species_names.append(file.split("-")[0])#Pugnoides weeksi fig28-3-c=Early Permian=Girty-1927.png
            img_lists.append(os.path.join(root,file))
            total_views.append(file.split("=")[0].split('-')[-1])
            genus.append(file.split(" ")[0])
            img_names.append(file)
    return img_lists, species_names, total_views, img_names, genus

if __name__ == '__main__': 
    num_upper_semilandmarks = 15
    num_lower_semilandmarks = 15
    num_front_semilandmarks = 15
    num_lateral_semilandmarks = 20
    a_view = np.array([])# choose which view of images need to be analysed
    #path = r'E:\Python_learning\PyTorch\U-2-Net-master\test_data\u2net_results_Ammonoidea'#########################
    path = r'E:\Python_learning\PyTorch\U-2-Net-master\test_data\u2net_ammonoid'
    img_lists, species_names, total_views, img_names, genus = get_list(path)
    for path in img_lists:
        coordinates = get_coordinates(path)
        four_landmarks_points = four_landmarks(coordinates)
        alandmarks = dorsal_view(four_landmarks_points,coordinates,num_upper_semilandmarks,num_lower_semilandmarks)
        a_view = np.append(a_view,alandmarks.reshape(alandmarks.shape[0] * 2))
    print("extraction finished")
    workbook = xlwt.Workbook(encoding = 'utf-8')
    worksheet = workbook.add_sheet("Ammonoidea")
    worksheet.write(0,0,label='img_name')
    worksheet.write(0,1,label='taxon')
    worksheet.write(0,2,label='genus')
    worksheet.write(0,3,label='view')
    
    for i in range(len(img_lists)):
        worksheet.write(i+1,0,label=img_names[i][:-4])
        worksheet.write(i+1,1,label=species_names[i])
        worksheet.write(i+1,2,label=genus[i])
        worksheet.write(i+1,3,label=total_views[i])
    workbook.save("./landmarks/ammonoid/Ammonoidea_lists.xls")#img_lists, species_names, total_views, img_names, genus 
    np.savetxt('./landmarks/ammonoid/Ammonoidea_landmarks.csv',a_view.reshape(int(a_view.shape[0]/(num_lateral_semilandmarks*4+4)),(num_lateral_semilandmarks*4+4)), delimiter=",",fmt='%d')
