/*
 * Emily Branson egb17c@acu.edu
 * 10/27/2021
 * 
 * This program takes the FST survey data and ideal FST reference point positions
 * to first create alignment matrices that transfer the ideal positions to the
 * measured positions, and rotate the ideal axes to the measured axes.
 * 
 * These matrices are outputted as ASCII files
 * 
 * This program also outputs ASCII files of the ideal, measured, and transferred
 * positions for every reference point of each module
 * 
 * This program also plots histograms of the ideal, measured, and transferred positions
 * with red = ideal, blue = measured, and green = transferred
 * 
 * Any piece of the code which would create and save a file is commented out
 * 
 */



#include "TSystem.h"



//use these variables and findTheta to find the angle difference between real and ideal
//vectors (middle to corner)
double r_real, r_ideal, difx_real, dify_real, difx_ideal, dify_ideal, rdoti,
       mag_real, mag_ideal;

double findTheta (double x, double y, double ix, double iy, double mx, double my)
{
    //find the difference between the point and the middle in terms of x and y
    difx_real = x - mx;
    dify_real = y - my;
    difx_ideal = ix - mx;
    dify_ideal = iy - my;
    
    //find the magnitudes of the vectors using x and y difference
    mag_real = sqrt(pow(difx_real,2) + pow(dify_real,2));
    mag_ideal = sqrt(pow(difx_ideal,2) + pow(dify_ideal,2));
    
    //find the angle difference between the two vectors
    rdoti = (difx_real * difx_ideal) + (dify_real * dify_ideal);
    return acos(rdoti/(mag_real*mag_ideal));
}

//application of rotation matrix (after transfer) to a point
double rotateX(double cosine, double sine, double xpos, double ypos, double xdif, double ydif, double midx,double midy)
{
         
    return midx + (xpos + xdif - midx)*cosine - (ypos + ydif - midy)*sine;
            
}

double rotateY(double cosine, double sine, double xpos, double ypos, double xdif, double ydif, double midx,double midy)
{
         
    return midy + (xpos + xdif - midx)*sine + (ypos + ydif - midy)*cosine;
            
}


//main function
void allMatricesInner()
{
    //This part makes an array of all files in the "STAR_files" directory
    const char* inDir = "/home/user/Documents/STAR_files/Inner";

    char* dir = gSystem->ExpandPathName(inDir);
    void* dirp = gSystem->OpenDirectory(dir);

    const char* entry;
    const char* filename[100];
    Int_t n = 0;
    TString str;

    //be sure to exclude the ideal data file from this array
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        str = entry;
        if(str.EndsWith(".txt")  && !(strncmp(str, "STAR", strlen("STAR"))))
            filename[n++] = entry;
    }
    

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...// 
   
   //make a tuple for each row, each tuple stores X,Y,Z for the row
   
    TNtuple row1("row1","Origin Center Point of Pin ID", "X:Y:Z");
    TNtuple row2("row2","Rotation Center Point of Pin ID", "X:Y:Z");
    TNtuple row3("row3","Lower Left Bias Intersection", "X:Y:Z");
    TNtuple row4("row4","Lower Right Bias Intersection", "X:Y:Z");
    TNtuple row5("row5","Left Pad Array Intersection Point", "X:Y:Z");
    TNtuple row6("row6","Left Middle Pad Array Intersection Point", "X:Y:Z");
    TNtuple row7("row7","Right Middle Pad Array Intersection Point", "X:Y:Z");
    TNtuple row8("row8","Right Pad Array Intersection Point", "X:Y:Z");
    TNtuple row9("row9","Bottom 'F' Corner Intersection", "X:Y:Z");

    double positions[50][12][2]; //52 files, 12 rows, 2 values (x then y)
    
    //this part loops through the array of filenames and opens every file one at a time
    for (int i = 0; i < n; i++){
    
        fstream file;
        string line;
        file.open(filename[i], ios::in);//make sure all files are in the same directory as the macro
        
        double x, y, z;
        
        
        if(file.is_open()){
            int count = 1;
            //go through every line of the file, incrementing "count" each time
            while(getline(file,line)){
                //set x, y, and z values
                if(line.substr(0,1) == "|" && line.substr(1,1) == "X"){
                
                    if(!line.substr(8,9).empty() && !line.substr(26,9).empty() && !line.substr(44,9).empty()){
                        x = stod(line.substr(8,9));
                        y = stod(line.substr(26,9));
                        z = stod(line.substr(44,9));
                    }
                    
                    //add the x and y positions to the "positions" array. Use these values later
                    //to calculate matrices
                    positions[i][count][0] = x;
                    positions[i][count][1] = y;
                    
                //use this switch to put the values from the row into the matching tuple
                    switch(count){
                        case 1:
                            row1.Fill(x,y,z);
                            break;
                        case 2:
                            row2.Fill(x,y,z);
                            break;
                        case 3:
                            row3.Fill(x,y,z);
                            break;
                        case 4:
                            row4.Fill(x,y,z);
                            break;
                        case 5:
                            row5.Fill(x,y,z);
                            break;
                        case 6:
                            row6.Fill(x,y,z);
                            break;
                        case 7:
                            row7.Fill(x,y,z);
                            break;
                        case 8:
                            row8.Fill(x,y,z);
                            break;
                        case 9:
                            row9.Fill(x,y,z);
                            break;
                    }
                    count ++;
                }
                
            }
            file.close();
        }
    }
    

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//
//This section makes and fills histograms with the ideal values



//delete all the old histograms to avoid errors
    delete gROOT->FindObject("h1x");
    delete gROOT->FindObject("h1xreal");
    delete gROOT->FindObject("h2x");
    delete gROOT->FindObject("h2xreal");
    delete gROOT->FindObject("h3x");
    delete gROOT->FindObject("h3xreal");
    delete gROOT->FindObject("h4x");
    delete gROOT->FindObject("h4xreal");
    delete gROOT->FindObject("h5x");
    delete gROOT->FindObject("h5xreal");
    delete gROOT->FindObject("h6x");
    delete gROOT->FindObject("h6xreal");
    delete gROOT->FindObject("h7x");
    delete gROOT->FindObject("h7xreal");
    delete gROOT->FindObject("h8x");
    delete gROOT->FindObject("h8xreal");
    delete gROOT->FindObject("h9x");
    delete gROOT->FindObject("h9xreal");
    
    delete gROOT->FindObject("h1y");
    delete gROOT->FindObject("h1yreal");
    delete gROOT->FindObject("h2y");
    delete gROOT->FindObject("h2yreal");
    delete gROOT->FindObject("h3y");
    delete gROOT->FindObject("h3yreal");
    delete gROOT->FindObject("h4y");
    delete gROOT->FindObject("h4yreal");
    delete gROOT->FindObject("h5y");
    delete gROOT->FindObject("h5yreal");
    delete gROOT->FindObject("h6y");
    delete gROOT->FindObject("h6yreal");
    delete gROOT->FindObject("h7y");
    delete gROOT->FindObject("h7yreal");
    delete gROOT->FindObject("h8y");
    delete gROOT->FindObject("h8yreal");
    delete gROOT->FindObject("h9y");
    delete gROOT->FindObject("h9yreal");
    
    delete gROOT->FindObject("h1xshift");
    delete gROOT->FindObject("h2xshift");
    delete gROOT->FindObject("h3xshift");
    delete gROOT->FindObject("h4xshift");
    delete gROOT->FindObject("h5xshift");
    delete gROOT->FindObject("h6xshift");
    delete gROOT->FindObject("h7xshift");
    delete gROOT->FindObject("h8xshift");
    delete gROOT->FindObject("h9xshift");
    
    delete gROOT->FindObject("h1yshift");
    delete gROOT->FindObject("h2yshift");
    delete gROOT->FindObject("h3yshift");
    delete gROOT->FindObject("h4yshift");
    delete gROOT->FindObject("h5yshift");
    delete gROOT->FindObject("h6yshift");
    delete gROOT->FindObject("h7yshift");
    delete gROOT->FindObject("h8yshift");
    delete gROOT->FindObject("h9yshift");
    
    delete gROOT->FindObject("c");
    
    
    
    //make the new histograms for the ideal measurement for each row
    //x or y is the ideal x or y, xreal or yreal is the measured x or y
    
    TH1D *h1x = new TH1D ("h1x","row1x ideal histogram",50,-1.2,1);
    TH1D *h1xreal = new TH1D ("h1xreal","X",50,-1.2,1);
    TH1D *h2x = new TH1D ("h2x","row2x ideal histogram",50,94.3,95.8);
    TH1D *h2xreal = new TH1D ("h2xreal","X",50,94.3,95.8);
    TH1D *h3x = new TH1D ("h3x","row3x ideal histogram",50,33.8,35.5);
    TH1D *h3xreal = new TH1D ("h3xreal","X",50,33.8,35.5);
    TH1D *h4x = new TH1D ("h4x","row4x ideal histogram",50,60,61.5);
    TH1D *h4xreal = new TH1D ("h4xreal","X",50,60,61.5);
    TH1D *h5x = new TH1D ("h5x","row5x ideal histogram",50,14.5,16);
    TH1D *h5xreal = new TH1D ("h5xreal","X",50,14.5,16);
    TH1D *h6x = new TH1D ("h6x","row6x ideal histogram",50,36,37.5);
    TH1D *h6xreal = new TH1D ("h6xreal","X",50,36,37.5);
    TH1D *h7x = new TH1D ("h7x","row7x ideal histogram",50,57.5,59.5);
    TH1D *h7xreal = new TH1D ("h7xreal","X",50,57.5,59.5);
    TH1D *h8x = new TH1D ("h8x","row8x ideal histogram",50,79,80.5);
    TH1D *h8xreal = new TH1D ("h8xreal","X",50,79,80.5);
    TH1D *h9x = new TH1D ("h9x","row9x ideal histogram",50,47,48.4);
    TH1D *h9xreal = new TH1D ("h9xreal","X",50,47,48.4);


    TH1D *h1y = new TH1D ("h1y","row1y ideal histogram",50,-1,1);
    TH1D *h1yreal = new TH1D ("h1yreal","Y",50,-1,1);
    TH1D *h2y = new TH1D ("h2y","row2y ideal histogram",50,-1,1);
    TH1D *h2yreal = new TH1D ("h2yreal","Y",50,-1,1);
    TH1D *h3y = new TH1D ("h3y","row3y ideal histogram",50,-254,-252);
    TH1D *h3yreal = new TH1D ("h3yreal","Y",50,-254,-252);
    TH1D *h4y = new TH1D ("h4y","row4y ideal histogram",50,-253.4,-252.4);
    TH1D *h4yreal = new TH1D ("h4yreal","Y",50,-253.4,-252.4);
    TH1D *h5y = new TH1D ("h5y","row5y ideal histogram",50,-140,-138.2);
    TH1D *h5yreal = new TH1D ("h5yreal","Y",50,-140,-138.2);
    TH1D *h6y = new TH1D ("h6y","row6y ideal histogram",50,-137,-135.5);
    TH1D *h6yreal = new TH1D ("h6yreal","Y",50,-137,-135.5);
    TH1D *h7y = new TH1D ("h7y","row7y ideal histogram",50,-137,-135.5);
    TH1D *h7yreal = new TH1D ("h7yreal","Y",50,-137,-135.5);
    TH1D *h8y = new TH1D ("h8y","row8y ideal histogram",50,-140,-138);
    TH1D *h8yreal = new TH1D ("h8yreal","Y",50,-140,-138);
    TH1D *h9y = new TH1D ("h9y","row9y ideal histogram",50,-252.5,-250.5);
    TH1D *h9yreal = new TH1D ("h9yreal","Y",50,-252.5,-250.5);
    
    
    //create new transferred ("shifted") histograms which will be used later
    TH1D *h1xshift = new TH1D ("h1xshift","X",50,-1.2,1);
    TH1D *h2xshift = new TH1D ("h2xshift","X",50,94.3,95.8);
    TH1D *h3xshift = new TH1D ("h3xshift","X",60,33.8,35.5);
    TH1D *h4xshift = new TH1D ("h4xshift","X",65,60,61.5);
    TH1D *h5xshift = new TH1D ("h5xshift","X",50,14.5,16);
    TH1D *h6xshift = new TH1D ("h6xshift","X",50,36,37.5);
    TH1D *h7xshift = new TH1D ("h7xshift","X",60,57.5,59.5);
    TH1D *h8xshift = new TH1D ("h8xshift","X",50,79,80.5);
    TH1D *h9xshift = new TH1D ("h9xshift","X",55,47,48.4);

             
    TH1D *h1yshift = new TH1D ("h1yshift","Y",50,-1,1);
    TH1D *h2yshift = new TH1D ("h2yshift","Y",50,-1,1);
    TH1D *h3yshift = new TH1D ("h3yshift","Y",100,-254,-252);
    TH1D *h4yshift = new TH1D ("h4yshift","Y",75,-253.4,-252.4);
    TH1D *h5yshift = new TH1D ("h5yshift","Y",100,-140,-138.2);
    TH1D *h6yshift = new TH1D ("h6yshift","Y",100,-137,-135.5);
    TH1D *h7yshift = new TH1D ("h7yshift","Y",75,-137,-135.5);
    TH1D *h8yshift = new TH1D ("h8yshift","Y",60,-140,-138);
    TH1D *h9yshift = new TH1D ("h9yshift","Y",60,-252.5,-250.5);
    
 
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//     
    //this section gets the values for the ideal histograms
    
    //open the ideal file and declare variables
    fstream file;
    string line;
    file.open("ideal_survey_inner.txt", ios::in);
    
    double x, y;
    double idealVals[18];
    double iTLx, iTLy, iTRx, iTRy, iBLx, iBLy, iBRx, iBRy, iRotx, iRoty, iOriginx, iOriginy,
           iLMx, iLMy, iRMx, iRMy, iBFx, iBFy;
    //open the file with the ideal values and fill the histograms with the data
    if(file.is_open()){
        int count = 1;
        while(getline(file,line)){
            if(line.substr(0,1) == "|" && line.substr(1,1) == "X"){
                
                if(!line.substr(9,9).empty() && !line.substr(26,9).empty()){
                    x = stod(line.substr(9,9));
                    y = stod(line.substr(29,9));
                }
                
                //store the x and y values in the array with x in first 9 indices, y in last 9
                idealVals[count - 1] = x;
                idealVals[count + 8] = y;
            
                switch(count){
                    case 1:
                        h1x->Fill(x);
                        h1y->Fill(y);
                        iOriginx = x;
                        iOriginy = y;
                        break;
                    case 2:
                        h2x->Fill(x);
                        h2y->Fill(y);
                        iRotx = x;
                        iRoty = y;
                        break;
                    case 3:
                        h3x->Fill(x);
                        h3y->Fill(y);
                        iBLx = x;
                        iBLy = y;
                        break;
                    case 4:
                        h4x->Fill(x);
                        h4y->Fill(y);
                        iBRx = x;
                        iBRy = y;
                        break;
                    case 5:
                        h5x->Fill(x);
                        h5y->Fill(y);
                        iTLx = x;
                        iTLy = y;
                        break;
                    case 6:
                        h6x->Fill(x);
                        h6y->Fill(y);
                        iLMx = x;
                        iLMy = y;
                        break;
                    case 7:
                        h7x->Fill(x);
                        h7y->Fill(y);
                        iRMx = x;
                        iRMy = y;
                        break;
                    case 8:
                        h8x->Fill(x);
                        h8y->Fill(y);
                        iTRx = x;
                        iTRy = y;
                        break;
                    case 9:
                        h9x->Fill(x);
                        h9y->Fill(y);
                        iBFx = x;
                        iBFy = y;
                        break;
                }
                count ++;
            }
            
        }
        file.close();
    }

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//
//this section is for calculating the matrices

    
    //declare arrays of 50 values for each of the 50 modules
    //origin is origin, rot is rotation point, BL = bottom left, BR = bottom right
    //TL = top left, TR = top right
    double OriginX[50], Rotx[50], Roty[50], BLx[50], BLy[50], BRx[50], BRy[50], TLx[50], TLy[50], TRx[50], TRy[50];
    
    //assign the variable names to the correct values from the data
    for (int a = 0; a < 50; a++){
        OriginX[a] = positions[a][1][0];
        Rotx[a] = positions[a][2][0];
        Roty[a] = positions[a][2][1];
        BLx[a] = positions[a][3][0];
        BLy[a] = positions[a][3][1];
        BRx[a] = positions[a][4][0];
        BRy[a] = positions[a][4][1];
        TLx[a] = positions[a][5][0];
        TLy[a] = positions[a][5][1];
        TRx[a] = positions[a][8][0];
        TRy[a] = positions[a][8][1];
    }
    

    //Mx = middle x, My = middle y. Middle 1 is calculated from the top of the sensor
    //middle 2 is calculated from the bottom of the sensor. MAV is the average of these
    double Mx1[50], Mx2[50], My1[50], My2[50], MxAV[50], MyAV[50];
    
    //calculate ideal middle
    double iMx1 = (iTRx + iTLx)/2;
    double iMx2 = (iBRx + iBLx)/2;
    double iMy1 = (iTLy + iBLy)/2;
    double iMy2 = (iTRy + iBRy)/2;
    double iMxAV = (iMx1 + iMx2)/2;
    double iMyAV = (iMy1 + iMy2)/2;
    
    
    double xdif[50], ydif[50], yrotdif[50], xrotdif[50], thetarot[50];
   
    double thetaTR, thetaTL, thetaBR, thetaBL;
    
           
    //xdif and ydif are found from the difference of the real and ideal middles
    for (int b = 0; b < 50; b++){
        Mx1[b] = (TRx[b] + TLx[b])/2;
        Mx2[b] = (BRx[b] + BLx[b])/2;
        My1[b] = (TLy[b] + BLy[b])/2;
        My2[b] = (TRy[b] + BRy[b])/2;

        MxAV[b] = (Mx1[b] + Mx2[b])/ 2;
        MyAV[b] = (My1[b] + My2[b])/2;

        xdif[b] = MxAV[b] - iMxAV;
        ydif[b] = MyAV[b] - iMyAV;

        //find the angle difference between real and ideal vector from middle to a corner
        //do this for each corner
        thetaTL = findTheta(TLx[b], TLy[b], iTLx, iTLy, MxAV[b], MyAV[b]);
        thetaTR = findTheta(TRx[b], TRy[b], iTRx, iTRy, MxAV[b], MyAV[b]);
        thetaBL = findTheta(BLx[b], BLy[b], iBLx, iBLy, MxAV[b], MyAV[b]);
        thetaBR = findTheta(BRx[b], BRy[b], iBRx, iBRy, MxAV[b], MyAV[b]);
        
        //average the angle difference for each corner
        thetarot[b] = (thetaTL + thetaTR + thetaBL + thetaBR)/4;
        
    }
    

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//    
//this section outputs files with the matrices
    
    /*
    ofstream myfile;
    string matricesFileName;
    
    for (int c = 0; c < n; c++){
        
        matricesFileName = filename[c];
        matricesFileName = "matrices_" + matricesFileName;
        myfile.open(matricesFileName);
        
        myfile << filename[c] << "\n\n" << endl;
        
        myfile << "Transfer matrix: \n";
        myfile << "     |       0       |       1       |\n";
        myfile << "--------------------------------------\n";
        myfile << "   0 |   "<<xdif[c]<<"  |   "<<ydif[c]<<"\n";
        
        myfile << "\nRotation matrix: \n";
        myfile << "     |       0       |       1       |\n";
        myfile << "--------------------------------------\n";
        myfile << "   0 |   "<<cos(thetarot[c]/2)<<"    |   "<<sin(thetarot[c]/2)<<"\n";
        myfile << "   1 |   "<<-sin(thetarot[c]/2)<<" |   "<<cos(thetarot[c]/2)<<"\n";
        
        cout << "created file   " << matricesFileName << endl;
        myfile.close();
    }
    */
    
       
//this next section transfers the ideal to real, and plots the real, ideal, and transferred positions together
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//    

    double sine, cosine;
    double posValues[18][50];
    double h1xtfer, h2xtfer, h3xtfer, h4xtfer, h5xtfer, h6xtfer, h7xtfer, h8xtfer, h9xtfer;
    double h1ytfer, h2ytfer, h3ytfer, h4ytfer, h5ytfer, h6ytfer, h7ytfer, h8ytfer, h9ytfer;
    
    //this for loop applies the transfer and rotation matrices to the original histograms
    for (int p = 0; p < n; p++){
        
        cosine = cos(thetarot[p]/2);
        sine = sin(thetarot[p]/2);
        
        //histograms go in this order: iOrigin, iRot, iBL, iBR, iTL, iLM, iRM, iTR, iBF
        
        
        //rotate the ideal to the real
        
        //X values
        h1xtfer=rotateX(cosine, sine, iOriginx, iOriginy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[0][p] = h1xtfer;
        h1xshift->Fill(h1xtfer);
        
        h2xtfer = rotateX(cosine, sine, iRotx, iRoty, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[0][p] = h2xtfer;
        h2xshift->Fill(h2xtfer);
        
        h3xtfer = rotateX(cosine, sine, iBLx, iBLy, xdif[p], ydif[p], MxAV[p], MyAV[p]); 
        posValues[2][p] = h3xtfer;
        h3xshift->Fill(h3xtfer);
        
        h4xtfer = rotateX(cosine, sine, iBRx, iBRy, xdif[p], ydif[p], MxAV[p], MyAV[p]); 
        posValues[3][p] = h4xtfer;
        h4xshift->Fill(h4xtfer);
        
        h5xtfer = rotateX(cosine, sine, iTLx, iTLy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[4][p] = h5xtfer;
        h5xshift->Fill(h5xtfer);
        
        h6xtfer = rotateX(cosine, sine, iLMx, iLMy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[5][p] = h6xtfer;
        h6xshift->Fill(h6xtfer);
        
        h7xtfer = rotateX(cosine, sine, iRMx, iRMy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[6][p] = h7xtfer;
        h7xshift->Fill(h7xtfer);
        
        h8xtfer = rotateX(cosine, sine, iTRx, iTRy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[7][p] = h8xtfer;
        h8xshift->Fill(h8xtfer);
        
        h9xtfer = rotateX(cosine, sine, iBFx, iBFy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[8][p] = h9xtfer;
        h9xshift->Fill(h9xtfer);
        
        ////////////////////////////////////////////////////////////////////////////////
        
        //Y values
        h1ytfer = rotateY(cosine, sine, iOriginx, iOriginy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[9][p] = h1ytfer;
        h1yshift->Fill(h1ytfer);
        
        h2ytfer = rotateY(cosine, sine, iRotx, iRoty, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[10][p] = h2ytfer;
        h2yshift->Fill(h2ytfer);
        
        h3ytfer = rotateY(cosine, sine, iBLx, iBLy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[11][p] = h3ytfer;
        h3yshift->Fill(h3ytfer);
        
        h4ytfer = rotateY(cosine, sine, iBRx, iBRy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[12][p] = h4ytfer;
        h4yshift->Fill(h4ytfer);
        
        h5ytfer = rotateY(cosine, sine, iTLx, iTLy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[13][p] = h5ytfer;
        h5yshift->Fill(h5ytfer);
        
        h6ytfer = rotateY(cosine, sine, iLMx, iLMy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[14][p] = h6ytfer;
        h6yshift->Fill(h6ytfer);
        
        h7ytfer = rotateY(cosine, sine, iRMx, iRMy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[15][p] = h7ytfer;
        h7yshift->Fill(h7ytfer);
        
        h8ytfer = rotateY(cosine, sine, iTRx, iTRy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[16][p] = h8ytfer;
        h8yshift->Fill(h8ytfer);
        
        h9ytfer = rotateY(cosine, sine, iBFx, iBFy, xdif[p], ydif[p], MxAV[p], MyAV[p]);
        posValues[17][p] = h9ytfer;
        h9yshift->Fill(h9ytfer);
    }
    
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//
//this section outputs files which give the ideal measured, and transferred values
//for each reference point
    
    /*
    ofstream myfile;
    string positionsFileName;
    
    string rowNames[9] = {"Rotation Center Point of Pin ID                 ",
                     "Left Sensor Lower Left  Bias Intersection       ",
                     "Left Sensor Lower Right Bias Intersection       ",
                     "Right Sensor Lower Left Bias Intersection       ",
                     "Right Sensor Lower Right Bias Intersection      ",
                     "Left Sensor Left Pad Array Intersection Point   ",
                     "Left Sensor Right Pad Array Intersection Point  ",
                     "Right Sensor Left Pad Array Intersection Point  ",
                     "Right Sensor Right Pad Array Intersection Point "};
    
    for (int s = 0; s < n; s++){
        
        positionsFileName = filename[s];
        positionsFileName = "positions_" + positionsFileName;
        myfile.open(positionsFileName);
        
        myfile << filename[s] << "\n\n" << endl;
        
        myfile << "IDEAL                    MEASURED                TRANSFERRED" << endl;
        
        for (int t = 0; t < 9; t++){
            
            myfile << fixed << "(" << idealVals[t] << ", " << setprecision(4) << idealVals[t + 9] << setprecision(4) << ")" << setw(3) << "|" << setw(3);
            myfile << fixed << "(" << positions[s][t][0] << ", " << setprecision(4) << positions[s][t][1] << ")" << setprecision(4) << setw(3) << "|" << setw(3);
            myfile << fixed << "(" <<posValues[t][s] << ", " << setprecision(4) << posValues[t + 9][s] << setprecision(4) << ")" << setw(3) << "|" << setw(3);
            myfile << rowNames[t] << "\n" << endl;
            
        }
        myfile.close();
        
        
    }
    
    */
    
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//

    //scale the ideal values on the plots with real measurements, and make the ideal lines red
    //x histograms
    h1x->Scale(25);
    h1x->SetLineColor(2);
    h2x->Scale(5.5);
    h2x->SetLineColor(2);
    h3x->Scale(5);
    h3x->SetLineColor(2);
    h4x->Scale(4.5);
    h4x->SetLineColor(2);
    h5x->Scale(3.5);
    h5x->SetLineColor(2);
    h6x->Scale(3.5);
    h6x->SetLineColor(2);   
    h7x->Scale(4.5);
    h7x->SetLineColor(2);  
    h8x->Scale(3.5);
    h8x->SetLineColor(2);   
    h9x->Scale(4.5);
    h9x->SetLineColor(2);   
    
    //y histograms
    h1y->Scale(22.5);
    h1y->SetLineColor(2);
    h2y->Scale(18.5);
    h2y->SetLineColor(2);
    h3y->Scale(5.5);
    h3y->SetLineColor(2);
    h4y->Scale(3.5);
    h4y->SetLineColor(2);
    h5y->Scale(7.5);
    h5y->SetLineColor(2);
    h6y->Scale(4);
    h6y->SetLineColor(2);
    h7y->Scale(4.5);
    h7y->SetLineColor(2);
    h8y->Scale(5.5);
    h8y->SetLineColor(2);
    h9y->Scale(5);
    h9y->SetLineColor(2);

    
    //set the linecolor of the shifted plots to green
    h1xshift->SetLineColor(kGreen + 1);
    h1yshift->SetLineColor(kGreen + 1);    
    h2xshift->SetLineColor(kGreen + 1);
    h2yshift->SetLineColor(kGreen + 1);
    h3xshift->SetLineColor(kGreen + 1);
    h3yshift->SetLineColor(kGreen + 1);
    h4xshift->SetLineColor(kGreen + 1);
    h4yshift->SetLineColor(kGreen + 1);
    h5xshift->SetLineColor(kGreen + 1);
    h5yshift->SetLineColor(kGreen + 1);
    h6xshift->SetLineColor(kGreen + 1);
    h6yshift->SetLineColor(kGreen + 1);
    h7xshift->SetLineColor(kGreen + 1);
    h7yshift->SetLineColor(kGreen + 1);
    h8xshift->SetLineColor(kGreen + 1);
    h8yshift->SetLineColor(kGreen + 1);
    h9xshift->SetLineColor(kGreen + 1);
    h9yshift->SetLineColor(kGreen + 1);
    
    


//this section plots the ideal, measured, and transferred values for every module
//ideal is red, measured is blue, transferred is green
  
    
    TCanvas *c2 = new TCanvas;
    
    for(int k = 1; k <= 11; k++){
        switch(k){
            case 1:
                row1.Draw("X>>h1xreal");
                h1x->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Origin_Center_Point_of_Pin_ID.png");
                h1xshift->Draw("SAME");
                
                row1.Draw("Y>>h1yreal");
                h1y->Draw("SAME");
                h1yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Origin_Center_Point_of_Pin_ID.png");
                break;
            case 2:
                row2.Draw("X>>h2xreal");
                h2x->Draw("SAME");
                h2xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Rotation_Center_Point_of_Pin_ID.png");
                
                row2.Draw("Y>>h2yreal");
                h2y->Draw("SAME");
                h2yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Rotation_Center_Point_of_Pin_ID.png");
                break;
            case 3:
                row3.Draw("X>>h3xreal");
                h3x->Draw("SAME");
                h3xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Lower_Left_Bias_Intersection.png");
                
                row3.Draw("Y>>h3yreal");
                h3y->Draw("SAME");
                h3yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Lower_Left_Bias_Intersection.png");
                break;
            case 4:
                row4.Draw("X>>h4xreal");
                h4x->Draw("SAME");
                h4xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Lower_Right_Bias_Intersection.png");
                
                row4.Draw("Y>>h4yreal");
                h4y->Draw("SAME");
                h4yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Lower_Right_Bias_Intersection.png");
                break;
            case 5:
                row5.Draw("X>>h5xreal");
                h5x->Draw("SAME");
                h5xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Left_Pad_Array_Intersection_Point.png");
                
                row5.Draw("Y>>h5yreal");
                h5y->Draw("SAME");
                h5yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Left_Pad_Array_Intersection_Point.png");
                break;
            case 6:
                row6.Draw("X>>h6xreal");
                h6x->Draw("SAME");
                h6xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Left_Middle_Pad_Array_Intersection_Point.png");
                
                row6.Draw("Y>>h6yreal");
                h6y->Draw("SAME");
                h6yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Left_Middle_Pad_Array_Intersection_Point.png");
                break;
            case 7:
                row7.Draw("X>>h7xreal");
                h7x->Draw("SAME");
                h7xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Right_Middle_Pad_Array_Intersection_Point.png");
                
                row7.Draw("Y>>h7yreal");
                h7y->Draw("SAME");
                h7yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Right_Middle_Pad_Array_Intersection_Point.png");
                break;
            case 8:
                row8.Draw("X>>h8xreal");
                h8x->Draw("SAME");
                h8xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Right_Pad_Array_Intersection_Point.png");
                
                row8.Draw("Y>>h8yreal");
                h8y->Draw("SAME");
                h8yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Right_Pad_Array_Intersection_Point.png");
                break;
            case 9:
                row9.Draw("X>>h9xreal");
                h9x->Draw("SAME");
                h9xshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/X_Bottom_'F'_Corner_Intersection.png");
                
                row9.Draw("Y>>h9yreal");
                h9y->Draw("SAME");
                h9yshift->Draw("SAME");
//                 c2->SaveAs("/home/user/Documents/STAR_files/Inner/innerTransfer/Y_Bottom_'F'_Corner_Intersection.png");
                break;
            }
    }
    
    //delete c2;
}


