# ParsimoniousKClub
This code accompanies the paper "Parsimonious formulations for low-diameter clusters" and is written in C++. If you want to use or cite this code, please cite the paper:   

            @misc{buchanan2018,
            title={Parsimonious formulations for low-diameter clusters},
            author={Buchanan, Austin and Salemi, Hosseinali},
            year={2018},
            note={Preprint available at \url{http://www.optimization-online.org/DB_HTML/2017/09/6196.html}}
            } 


Compiling thie code
--------------------
#### How to run the code for "Real-Life Instances (DIMACS-10)"

1.0. Download all .cpp and .h files. 
2.0. Start a new C++ project and name it "Kclub". Add all .cpp files to "Kclub" Source files and all .h files to its Header.   

3.0. Build the solution (CTRL+SHIFT+B). 

4.0. Open Windows Command Prompt.

5.0. Find "Release" folder in the "Kclub" folder. Type " cd [folder address of "Release"] " and press Enter. 

6.1. To see the results for "Recursive (R) formulation", type " Kclub.exe Veremyev dimacs [folder address of the instance] \ [name of instance].graph [value of k] " and press Enter. 

6.2. To see the results for "Canonical hypercube cut (CHC) formulation", type " Kclub.exe CHC dimacs [folder address of the instance] \ [name of instance].graph [value of k] " and press Enter.

6.3. To see the results for "Cut-like formulation (CUT)", type " Kclub.exe CutLike dimacs [folder address of the instance] \ [name of instance].graph [value of k] " and press Enter.
 
6.4. To see the results for "Path-like formulation (PATH)" where k=3 type " Kclub.exe 3ClubPathLike dimacs [folder address of the instance] \ [name of instance].graph " and press Enter.

6.5. To see the results for "Path-like formulation (PATH)" where k=4 type " Kclub.exe 4ClubPathLike dimacs [folder address of the instance] \ [name of instance].graph " and press Enter.

6.6. To see the results for "CN formulation" where k=2, type " Kclub.exe S2 dimacs [folder address of the instance] \ [name of instance].graph" and press Enter.

6.7. To see the results for "CN+ formulation" where k=2, type " Kclub.exe CutLike dimacs [folder address of the instance] \ [name of instance].graph 2" and press Enter.


#### How to run the code for "Synthetic instances "

1.0. Download all .cpp and .h files.

2.0. Start a new C++ project and name it "Kclub". Add all .cpp files to "Kclub" Source files and all .h files to its Header.   

3.0. Build the solution (CTRL+SHIFT+B). 

4.0. Open Windows Command Prompt.

5.0. Find "Release" folder in the "Kclub" folder. Type " cd [folder address of "Release"] " and press Enter. 

6.1. To see the results for "Recursive (R) formulation", type " Kclub.exe Veremyev snap_d [folder address of the instance] \ [name of instance].txt [value of k] " and press Enter. 

6.2. To see the results for "Canonical hypercube cut (CHC) formulation", type " Kclub.exe CHC snap_d [folder address of the instance] \ [name of instance].txt [value of k] " and press Enter. 

6.3. To see the results for "Cut-like formulation", type " Kclub.exe CutLike snap_d [folder address of the instance] \ [name of instance].txt [value of k] " and press Enter. 

6.4. To see the results for "Path-like formulation" where k=3, type " Kclub.exe 3ClubPathLike snap_d [folder address of the instance] \ [name of instance].txt " and press Enter. 

6.5. To see the results for "Path-like formulation" where k=4, type " Kclub.exe 4ClubPathLike snap_d [folder address of the instance] \ [name of instance].txt " and press Enter. 



Terms and conditions
--------------------
Copyright (c) 2018 Austin L. Buchanan, Hosseinali Salemi. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
