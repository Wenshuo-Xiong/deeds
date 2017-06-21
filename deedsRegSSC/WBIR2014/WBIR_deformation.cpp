/*
 Copyright (c) 2014, Mattias P. Heinrich
 Contact: heinrich(at)imi.uni-luebeck.de
 www.mpheinrich.de
 
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

/*
 If you use this implementation or parts of it please cite:
 
 “Non-parametric Discrete Registration with Convex Optimisation.”
 by Mattias P. Heinrich, Bartlomiej W. Papież, Julia A. Schnabel, Heinz Handels
 Biomedical Image Registration - WBIR 2014, LNCS 8454, Springer, pp 51-61
 
 and
 
 "Towards Realtime Multimodal Fusion for Image- Guided Interventions using Self-Similarities"
 by Mattias P. Heinrich, Mark Jenkinson, Bartlomiej W. Papiez, Sir Michael Brady, and Julia A. Schnabel
 Medical Image Computing and Computer-Assisted Intervention - MICCAI 2013, LNCS 8149, Springer
 
 
 Tested with g++ on Mac OS X and Linux Ubuntu, compile with:
 
 g++ WBIR_deformation.cpp -O3 -lpthread -msse4.2 -o defWBIR
 
 replace msse4.2 by your current SSE version if needed
 for Windows you might need MinGW or CygWin

*/
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sys/time.h>
#include <math.h>
#include <inttypes.h>
#include <pthread.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

using namespace std;


#include "deeds2deformation.h"

int main(int argc, const char * argv[])
{
    
    if(argc<3||argv[1][1]=='h'){
        cout<<"==========================================================\n";
        cout<<"Usage (required input arguments):\n";
        cout<<"./defWBIR -F fixed.nii -O output\n";
        cout<<"==========================================================\n";
        return 1;
    }
    
    
    typedef pair<char,int> val;
    map<char,int> argin;
    argin.insert(val('F',0));
    argin.insert(val('O',1));

    // parsing the input
    int requiredArgs=0;
     char* fixedfile=new char[200];
     char* outputstem=new char[200];
    
    
    for(int k=1;k<argc;k++){
        if(argv[k][0]=='-'){
            if(argin.find(argv[k][1])==argin.end()){
                cout<<"Invalid option: "<<argv[k]<<" use -h for help\n";
            }
            switch(argin[argv[k][1]]){
                case 0:
                    sprintf(fixedfile,"%s",argv[k+1]);
                    requiredArgs++;
                    break;
                case 1:
                    sprintf(outputstem,"%s",argv[k+1]);
                    requiredArgs++;
                    break;
                default:
                    cout<<"Invalid option: "<<argv[k]<<" use -h for help\n";
                    break;
               }
        }
    }
    if(requiredArgs!=2){
        cout<<"Missing argmuents, use -h for help.\n";
    }
    
    
    deeds2def(fixedfile,outputstem);			 //只有一个fixed 和output

    
    return 0;
}

