This is a matlab implementation of the method described in:

"Minimum Barrier Salient Object Detection at 80 FPS", Jianming Zhang, 
Stan Sclaroff, Zhe Lin, Xiaohui Shen, Brian Price, Radomir Mech, ICCV, 2015

Contact: jimmie33@gmail.com

Prerequisite: OpenCV 2.4+

Usage:

1. Go to the folder "mex"
2. modify the opencv include and lib paths in "compile.m/compile_win.m" 
   (for Linux/Windows)
3. run "compile/compile_win" in matlab (for Linux/Windows)
4. Go to the root folder
5. run "demo"


This matlab implementation is provided for research purpose only. For fully
reproducing the results in our ICCV paper, please use the original Windows
executable program. 

The matlab implementation is slower than the window executable, mainly due 
to the morphological post-processing step. We use the highly optimized IPP
library for the morphological operations in our C++ implementation, which 
are much faster than the corresponding Matlab functions.  

