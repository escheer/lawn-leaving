These MATLAB functions are for constructing a dictionary of behavioural motifs and using the dictionary to make a phenotypic similarity network that can be displayed using, for example, Cytoscape.  If you use this code in a paper, please cite our paper describing this work:

André E.X. Brown, Eviatar I. Yemini, Laura Grundy, Tadas Jucikas, William R. Schafer (2013) A dictionary of behavioral motifs reveals clusters of genes affecting Caenorhabditis elegans locomotion. Proc. Nat. Acad. Sci. USA 110:791.  http://dx.doi.org/10.1073/pnas.1211447110

Everything you need to run the motifAnalysis_DEMO.m script should be included.  You may need to recompile the mex functions for your system.  The source code is included.  See http://www.mathworks.co.uk/support/tech-notes/1600/1605.html for instructions.

I am grateful for several contributed functions including the motif finding functions from Abdullah Mueen and Eamonn Keogh, timeseries2symbol from Eamonn Keogh, Jessica Lin, Stefano Lonardi, and Pranav Patel, the rdir function from Gus Brown, covshrinkKPM from Kevin Murphy, and the mRMR package from Hanchuan Peng.  These functions all have their own licenses included in their headers or, in the case of mRMR, in its readme file.

If you find any bugs or have any questions, please e-mail me at abrown@mrc-lmb.cam.ac.uk or aexbrown@gmail.com.


The MIT License

Copyright (c)  Medical Research Council 2013

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.