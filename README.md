# DeepMorph
Morphological analysis of the latest Permian to earliest Triassic marine fossils. 1.) We collected the fossil images from publications; 2) U2-Net is used to extract the binary images; 3) ammonoids_landmarks.py can extract the contours of fossils and we used geometric morphometrics to get landmarks and semi-landmarks, and coordinates will be saved for morphological analyze; 4) alignment of specimens is based on R package, i.e., geomorph, then we used disparity for morphological quantification. 

![image](https://github.com/XiaokangLiuCUG/DeepMorph/blob/main/Figure%201%20Schematic%20of%20pipeline.png)


## How to use DeepMorph
1. Prepare the DeepMorph environment based on the requirements.txt. The DeepMorph was implemented in Python (> 3.6.5), Pytorch (v. 1.10.2), and torchvision (v. 0.11.3+cu113).
 
3. U2-Net was  used to segment fossil areas by running the u2net_test.py. Please add your model and image path accordingly.
 ```
python .../u2net_test.py
```
3. Then, the contours of the fossil can be extracted using ammonoids_landmarks.py. ps change to your path of the binary image file obtained from the second step. Fossil contours will be saved as coordinates in a CSV file. Among them, four landmarks and 60 will be extracted as shown in the figure above (trait matrix). Here, I only used one standard view for each fossil specimen.  
 ```
python .../u2net_test.py 
```
4. Generalized Procrustes analysis (GPA) was conducted to align shells and remove the size effect, then we can calculate the disparity of different subsets, as demonstrated in PT_ammonoids_disparity.R.

## Citation
Liu X, Song H, Chu D, Dai, X., Wang, F., and Silvestro D. Heterogeneous selectivity and morphological evolution of marine clades during the Permian–Triassic mass extinction. Nature Ecology & Evolution, 2024: 8: 1248–1258
