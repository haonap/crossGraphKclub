# Code for computing cross-graph k-clubs
This repository contains C++ code used for computing cross-graph k-clubs in the article "Finding Conserved Low-Diameter Subgraphs in Social and Biological Networks" which has been submitted to Networks. If you wish to use or cite this code, please cite:
        
        @article{Hao2024Finding,
                author = {Hao Pan and Yajun Lu and Balabhaskar Balasundaram and Juan S. Borrero},
                journal = {Networks},
                month = {February},
                note = {Under Review},
                title = {Finding Conserved Low-Diameter Subgraphs in Social and Biological Networks},
                year = {2024}
        }

# Understanding and using the code
The code should be straightforward if you start to read from file main.cpp. Necessary comments have been added in the code for easiness of understanding. Descriptions are added at the top of each function in functions.cpp and classes.cpp.

In file main.cpp, the main function starts by reading in input parameters from file tasksCrossGraph.txt. This file comprises six entries: the number of graph collections (p), the index of the first graph, the index of the last graph, the instance name, k, and the method (which can either be CCF or PPCF). For instance, an entry might look like this: "2 1 2 graph_200_0.1 2 PPCF".

# Compilation and execution in Linux environment
1. Download or clone the repository to your machine. 
2. From terminal, go to the repository. 
3. Type "make" and hit enter to compile. 
4. Open tasksCrossGraph.txt file to configure parameters. 
5. In terminal, type "./main" and hit enter to execute. 


# MIT License

Copyright (c) 2024 Hao Pan, Yajun Lu, Balabhaskar Balasundaram, and Juan S. Borrero

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
