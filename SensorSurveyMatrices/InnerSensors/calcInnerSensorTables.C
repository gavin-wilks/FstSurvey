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

/* Gavin Wilks gwilks3@uic.edu
 * 02/04/2022
 *
 * Modifying Emily Branson's original program to properly calculate the misalignment
 * tables for the FST sensors in the frame of the half disk.
 *
 */


#include "TSystem.h"

//use these variables and findTheta to find the angle difference between real and ideal
//vectors (middle to corner)
double r_real, r_ideal, difx_real, dify_real, difx_ideal, dify_ideal, rdoti,
       mag_real, mag_ideal;

int a_moduleRot[36] = {1,0,1,0,1,0,1,0,1,0,1,0,
                       0,1,0,1,0,1,0,1,0,1,0,1,
                       1,0,1,0,1,0,1,0,1,0,1,0};

double lengthMidPin = TMath::Sqrt(TMath::Sqrt( 303.33000*303.33000 + 31.88000*31.88000 )*TMath::Sqrt( 303.33000*303.33000 + 31.88000*31.88000 ) - 47.71000*47.71000);

void rotatePoint (double &x, double &y, double angle)
{
  //Rotate a 2D point by an angle

  double xtemp = x;
  double ytemp = y;
  
  x = xtemp*TMath::Cos(angle) - ytemp*TMath::Sin(angle);
  y = xtemp*TMath::Sin(angle) + ytemp*TMath::Cos(angle);
}

double findCoordinateRotation (double rotmx, double rotmy)
{
  //Find the coorindate rotation after shifting rotation point into the frame where originpin = (0,0)

  double length = TMath::Sqrt( rotmx*rotmx + rotmy*rotmy );
  double denomi = ( rotmx*rotmx/rotmy + rotmy);

  return TMath::ASin( -length/denomi );
}

double findTheta (double ix, double iy, double mxpp, double mypp, double refpointx, double refpointy)
{
    //Find the angle of rotation to place sensor in correct orientation of FST Half
    double idelx = refpointx - ix;
    double idely = refpointy - iy;
    double mdelx = refpointx - mxpp;
    double mdely = refpointy - mypp;

    double ilength = TMath::Sqrt( idelx*idelx + idely*idely ); //length from alignment pin to ideal middle
    double mlength = TMath::Sqrt( mdelx*mdelx + mdely*mdely ); //length from alignment pin to measured middle

    double numer = (ilength/mlength) * (mdely - idely*mdelx/idelx);
    double denom = idelx + idely*idely/idelx;

    return TMath::ASin( numer/denom );
}

std::string a_moduleId[36] = {"FST-73","FST-75","FST-62","FST-72","FST-76","FST-59",
                              "FST-42","FST-54","FST-39","FST-63","FST-47","FST-28",
                              "FST-69","FST-78","FST-50","FST-27","FST-31","FST-25",
                              "FST-34","FST-30","FST-24","FST-33","FST-22","FST-51",
                              "FST-41","FST-65","FST-37","FST-67","FST-45","FST-64",
                              "FST-21","FST-48","FST-17","FST-43","FST-55","FST-61"};


//main function
void calcInnerSensorTables(std::string inputlist = "../list/InnerSensors.list")
{

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

    double positions[36][12][2]; //36 files, 12 rows, 2 values (x then y)

    int i_file = 0;
    if (!inputlist.empty())   // if input file is ok
    {
      std::cout << "Open input file list: " << inputlist.c_str() << std::endl;
      std::ifstream in(inputlist.c_str());  // input stream
      if(in)
      {
        std::cout << "input database file list is ok" << std::endl;
        char str[255];       // char array for each file name
        while(in)
        {
          in.getline(str,255);  // take the lines of the file list
          if(str[0] != 0)
          {
            std::string inputfile;
            inputfile = str;
            std::cout << "open file: " << inputfile.c_str() << std::endl;
            std::fstream file;
            std::string line;
            file.open(inputfile, ios::in);

            double x, y, z;
        
        
            if(file.is_open())
            {
              int count = 0;
              //go through every line of the file, incrementing "count" each time
              while(getline(file,line))
              {
                //set x, y, and z values
                if(line.substr(0,1) == "|" && line.substr(1,1) == "X")
                { 
                  if(!line.substr(8,9).empty() && !line.substr(26,9).empty() && !line.substr(44,9).empty())
                  {
                    x = stod(line.substr(8,9)); 
                    y = stod(line.substr(26,9));
                    z = stod(line.substr(44,9));
                  }
                  
                  //add the x and y positions to the "positions" array. Use these values later
                  //to calculate matrices
                  positions[i_file][count][0] = x;
                  positions[i_file][count][1] = y;
                        
                  //use this switch to put the values from the row into the matching tuple
                  switch(count){
                    case 0:
                      row1.Fill(x,y,z);
                      break;
                    case 1:
                      row2.Fill(x,y,z);
                      break;
                    case 2:
                      row3.Fill(x,y,z);
                      break;
                    case 3:
                      row4.Fill(x,y,z);
                      break;
                    case 4:
                      row5.Fill(x,y,z);
                      break;
                    case 5:
                      row6.Fill(x,y,z);
                      break;
                    case 6:
                      row7.Fill(x,y,z);
                      break;
                    case 7:
                      row8.Fill(x,y,z);
                      break;
                    case 8:
                      row9.Fill(x,y,z);
                      break;
                  }
                  count ++;
                }
                    
              }
              file.close();
            }
            i_file++;
          }
        }
      }
    }

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//
//This section makes and fills histograms with the ideal values

    
    //make the new histograms for the ideal measurement for each row
    //x or y is the ideal x or y, xreal or yreal is the measured x or y
    TGraph *gRotVLength = new TGraph();      
    TGraph *gOriginVLength = new TGraph();   

    TH2D *hist2D_before[9];
    hist2D_before[0] = new TH2D("h_mOriginDiffB","M-AI Origin Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[1] = new TH2D("h_mRotDiffB","M-AI Rot Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[2] = new TH2D("h_mBLDiffB","M-AI BL Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[3] = new TH2D("h_mBRDiffB","M-AI BR Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[4] = new TH2D("h_mTLDiffB","M-AI TL Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[5] = new TH2D("h_mLMDiffB","M-AI LM Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[6] = new TH2D("h_mRMDiffB","M-AI RM Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[7] = new TH2D("h_mTRDiffB","M-AI TR Before",50,-1.0,1.0,50,-1.0,1.0);
    hist2D_before[8] = new TH2D("h_mBFDiffB","M-AI BF Before",50,-1.0,1.0,50,-1.0,1.0);

    TH2D *hist2D[9];
    hist2D[0] = new TH2D("h_mOriginDiff","M-AI Origin",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[1] = new TH2D("h_mRotDiff","M-AI Rot",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[2] = new TH2D("h_mBLDiff","M-AI BL",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[3] = new TH2D("h_mBRDiff","M-AI BR",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[4] = new TH2D("h_mTLDiff","M-AI TL",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[5] = new TH2D("h_mLMDiff","M-AI LM",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[6] = new TH2D("h_mRMDiff","M-AI RM",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[7] = new TH2D("h_mTRDiff","M-AI TR",50,-1.0,1.0,50,-1.0,1.0);
    hist2D[8] = new TH2D("h_mBFDiff","M-AI BF",50,-1.0,1.0,50,-1.0,1.0);
 
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//     
    //this section gets the values for the ideal histograms
    
    //open the ideal file and declare variables
    std::fstream file;
    std::string line;
    file.open("../data/InnerSensors/ideal_survey_inner.txt", ios::in);
    
    double x, y;
    double idealVals[9][2];
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
                idealVals[count - 1][0] = x;
                idealVals[count - 1][1] = y;
            
                switch(count){
                    case 1:
                        iOriginx = x;
                        iOriginy = y;
                        break;
                    case 2:
                        iRotx = x;
                        iRoty = y;
                        break;
                    case 3:
                        iBLx = x;
                        iBLy = y;
                        break;
                    case 4:
                        iBRx = x;
                        iBRy = y;
                        break;
                    case 5:
                        iTLx = x;
                        iTLy = y;
                        break;
                    case 6:
                        iLMx = x;
                        iLMy = y;
                        break;
                    case 7:
                        iRMx = x;
                        iRMy = y;
                        break;
                    case 8:
                        iTRx = x;
                        iTRy = y;
                        break;
                    case 9:
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

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//

    std::fstream infile;
    infile.open("../data/InnerSensors/InnerSensorAlignmentPins.csv", ios::in);

    float RefPointX[36], RefPointY[36];
    std::string templine;
    float xtemp, ytemp;
    char ch;
    int scount = 0;
   
    if(infile.is_open())
    {
      while(std::getline(infile, templine))
      {
        std::stringstream str(templine);
        
        if(scount == 0)
        {
          scount++;
          std::cout << str.str() << endl;
          continue;
        }
  
        str >> xtemp; RefPointX[scount-1] = xtemp; str >> ch; 
        str >> ytemp; RefPointY[scount-1] = ytemp;                      
  
        std::cout << scount++ << endl; 
      }
    }

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//

//this section is for calculating the matrices

    
    //declare arrays of 50 values for each of the 50 modules
    //origin is origin, rot is rotation point, BL = bottom left, BR = bottom right
    //TL = top left, TR = top right
    double OriginX[36], OriginY[36], Rotx[36], Roty[36], BLx[36], BLy[36], BRx[36], BRy[36], TLx[36], TLy[36], TRx[36], TRy[36];
    
    //assign the variable names to the correct values from the data
    for (int a = 0; a < 36; a++){
        OriginX[a] = positions[a][0][0];
        OriginY[a] = positions[a][0][1];
        Rotx[a] = positions[a][1][0];
        Roty[a] = positions[a][1][1];
        BLx[a] = positions[a][2][0];
        BLy[a] = positions[a][2][1];
        BRx[a] = positions[a][3][0];
        BRy[a] = positions[a][3][1];
        TLx[a] = positions[a][4][0];
        TLy[a] = positions[a][4][1];
        TRx[a] = positions[a][7][0];
        TRy[a] = positions[a][7][1];
    }
    

    //Mx = middle x, My = middle y. Middle 1 is calculated from the top of the sensor
    //middle 2 is calculated from the bottom of the sensor. MAV is the average of these
    double Mx1[36], Mx2[36], My1[36], My2[36], MxAV[36], MyAV[36];

    
    //calculate ideal middle
    double iMx1 = (iTRx + iTLx)/2;
    double iMx2 = (iBRx + iBLx)/2;
    double iMy1 = (iTLy + iBLy)/2;
    double iMy2 = (iTRy + iBRy)/2;
    double iMxAV = (iMx1 + iMx2)/2;
    double iMyAV = (iMy1 + iMy2)/2;
    
    
    double xdif[36], ydif[36], thetarot[36];          
    double deltax[36], deltay[36], angle[36], moduleAngle[36], lengthOMDiff[36];// QA Variables;

    //xdif and ydif are found from the difference of the real and ideal middles
    for (int b = 0; b < 36; b++)
    {
      Mx1[b] = (TRx[b] + TLx[b])/2;
      Mx2[b] = (BRx[b] + BLx[b])/2;
      My1[b] = (TLy[b] + BLy[b])/2;
      My2[b] = (TRy[b] + BRy[b])/2;

      MxAV[b] = (Mx1[b] + Mx2[b])/2;
      MyAV[b] = (My1[b] + My2[b])/2;
     
      //////////////////////////////////////////////////////////////////////////////////////////////////
      deltax[b] = -OriginX[b]; //Will translate all points to the coordinate frame where originpin = (0,0)
      deltay[b] = -OriginY[b];
 
      std::cout << "Delta x = " << deltax[b] << " Delta y = " << deltay[b] << std::endl;
      std::cout << "Before translation RotPinx = " << Rotx[b] << " Rotpiny = " << Roty[b] << std::endl;
      std::cout << "Before translation middlex = " << MxAV[b] << " middley = " << MyAV[b] << std::endl; 
      double rotmxp = Rotx[b] + deltax[b]; //Transfer rotation point into the proper frame
      double rotmyp = Roty[b] + deltay[b];
      double mxp = MxAV[b] + deltax[b];       //Transfer sensor middle point to proper frame
      double myp = MyAV[b] + deltay[b];
        
      std::cout << "Ideal      : Middle x = " << std::setprecision(10) << iMxAV  << " Middle y = " << iMyAV << std::endl;
      std::cout << "After shift: Middle x = " << mxp << " Middle y = " << myp << std::endl;

      angle[b] = findCoordinateRotation(rotmxp, rotmyp); //Angle of rotation s.t. rotpin = (x'',0)
     
      std::cout << "Rotation angle for rotation pin = " << angle[b] << std::endl;
      std::cout << "Before rotation RotPinx = " << rotmxp << " Rotpiny = " << rotmyp << std::endl;

      rotatePoint(rotmxp, rotmyp, angle[b]); //Rotates rotpin to proper position
      rotatePoint(mxp, myp, angle[b]);       //Rotates sensor middle to proper position

      std::cout << "After rotation RotPinx = " << rotmxp << " Rotpiny = " << rotmyp << std::endl;
      //////////////////////////////////////////////////////////////////////////////////////////////////
 
      //////////////////////////////////////////////////////////////////////////////////////////////////
      //Shift points s.t. the ideal midpoint would be at (0, lengthMidPin) 
      double mxpp;
      double imxpp = 0.0;
      double ox, oy; //new 
      if(a_moduleRot[b] == 0) 
      {
        mxpp = -(mxp - iMxAV);
        ox = 47.71; //new
        rotmxp = -(rotmxp - idealVals[1][0]);
      }
      else
      {
        mxpp = mxp - iMxAV; 
        ox = -47.71; //new
        rotmxp = rotmxp - idealVals[1][0];
      }
      oy = lengthMidPin; //new
      double mypp = lengthMidPin + myp;
      double imypp = lengthMidPin + iMyAV;
      rotmyp = lengthMidPin;

      std::cout << "Place module on x = 0,  x = " << mxpp << ", y = " << mypp << std::endl;

      //Rotate midpoint to proper location
      moduleAngle[b] = -TMath::Pi()/12.0 - (b%12)*TMath::Pi()/6.0;
      std::cout << "Module angle: " << TMath::RadToDeg()*moduleAngle[b] << std::endl;
      rotatePoint(ox, oy, moduleAngle[b]); //new 
      rotatePoint(mxpp, mypp, moduleAngle[b]);
      rotatePoint(imxpp, imypp, moduleAngle[b]); 
      rotatePoint(rotmxp,rotmyp, moduleAngle[b]);

      std::cout << "Ideal midpoint location in STAR Frame: x = " << imxpp << ", y = " << imypp << std::endl;
      std::cout << "Measured midpoint location in STAR Frame: x = " << mxpp << ", y = " << mypp << std::endl;
      //////////////////////////////////////////////////////////////////////////////////////////////////

      std::cout << "Ideal origin in STAR Frame: x = " << ox << ", y = " << oy << std::endl;
      //Find angle of rotation in FST Half coordinate frame
      //thetarot[b] = findTheta(imxpp, imypp, mxpp, mypp, RefPointX[b], RefPointY[b]); //old
     
      //******QA******
      double idelx = ox - imxpp;
      double idely = oy - imypp;
      double mdelx = ox - mxpp;
      double mdely = oy - mypp;

      double ilength = TMath::Sqrt( idelx*idelx + idely*idely ); //length from alignment pin to ideal middle
      double mlength = TMath::Sqrt( mdelx*mdelx + mdely*mdely ); //length from alignment pin to measured middle

      lengthOMDiff[b] = TMath::Abs(ilength - mlength);
      //******QA******

      thetarot[b] = findTheta(imxpp, imypp, mxpp, mypp, ox, oy); //new
      std::cout << "Angle = " << thetarot[b] << std::endl;
     
      //Rotate the ideal midpoint of the sensor to match measured orientation 
      rotatePoint(imxpp, imypp, thetarot[b]);
      std::cout << "Rotated ideal midpoint: x = " << imxpp << ", y = " << imypp << std::endl;

      //Calculte shift required to place the rotated ideal point onto the measured 
      //midpoint position w.r.t. the proper '' coordinate frame origin
      xdif[b] = mxpp - imxpp;
      ydif[b] = mypp - imypp; 
      std::cout << "xdif = " << xdif[b] << "   ydif = " << ydif[b] << std::endl << std::endl;
    }
    

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//    
    for(int i = 0; i < 9; i++) // loop over measurement points
    {
      for(int m = 0; m < 36; m++) // loop over modules
      {
        double length; //used for a plot at the end

        double x = positions[m][i][0];
        double y = positions[m][i][1];
        
        std::cout << "Measured positions:                        x=" << x << ", y=" << y << std::endl; 
      
        x += deltax[m]; 
        y += deltay[m];
        rotatePoint(x, y, angle[m]);
 
        if(i == 1) length = x;

        std::cout << "Recentered measured positions:             x=" << x << ", y=" << y << std::endl; 

        double ix, iy;
        if(a_moduleRot[m] == 0) 
        {
          x = -(x - 47.71); // 47.71 is the x center of the module in the sensor frame
          ix = -(idealVals[i][0] - 47.71);              
        }
        else
        {
          x = x - 47.71; 
          ix = idealVals[i][0] - 47.71;
        }
        y = lengthMidPin + y; 
        iy = lengthMidPin + idealVals[i][1];

        std::cout << "Recentered measured STAR positions:        x=" << x << ", y=" << y << std::endl; 
        std::cout << "Recentered ideal    STAR positions:        x=" << ix << ", y=" << iy << std::endl; 

        rotatePoint(x, y, moduleAngle[m]);   //rotate measured point to proper STAR coordinates
        rotatePoint(ix, iy, moduleAngle[m]); //rotate ideal point to proper STAR coordinates

        std::cout << "Recentered rotated measured STAR positions:x=" << x << ", y=" << y << std::endl; 
        std::cout << "Recentered rotated ideal   STAR positions: x=" << ix << ", y=" << iy << std::endl; 
       
        hist2D_before[i]->Fill(x - ix, y - iy);

        rotatePoint(ix, iy, thetarot[m]); //rotate ideal point using angle from alignment matrices

        std::cout << "Apply alignment rotation:                  x=" << ix << ", y=" << iy << std::endl; 

        ix += xdif[m]; //translate ideal point using translations from alignment matrices         
        iy += ydif[m];

        std::cout << "Apply alignment shift:                     x=" << ix << ", y=" << iy << std::endl;
        std::cout << std::endl;

        hist2D[i]->Fill(x - ix, y - iy); //Measured x,y - Aligned Ideal x,y in STAR coordinate frame 
       
        if(i == 0)
        {
          gOriginVLength->SetPoint(m,lengthOMDiff[m],TMath::Sqrt((x-ix)*(x-ix)+(y-iy)*(y-iy)));
        }
        if(i == 1)
        {
          gRotVLength->SetPoint(m,95.42-length,TMath::Sqrt((x-ix)*(x-ix)+(y-iy)*(y-iy)));
        }
        
      }
    }

//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//    
//this section outputs files with the matrices
    
    std::ofstream outfile;
    std::string outputname = "./output/InnerAlignmentValues.txt";
    outfile.open(outputname, ios::out);
  
    outfile << "ModuleID Rotation(rad) TranslationX(mm) TranslationY(mm)" << std::endl;

    for(int imod = 0; imod < 36; imod++)
    {
      outfile << a_moduleId[imod] << " " << thetarot[imod] << " " << xdif[imod]/10.0 << " " << ydif[imod]/10.0 << std::endl;
    }
 
    outfile.close();
          
//......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo......ooOO00OOoo...//
    TCanvas *c1 = new TCanvas("c1","c1",10,10,900,900);
    c1->Divide(3,3);
    for(int ip = 0; ip < 9; ip++)
    {
      c1->cd(ip+1)->SetLeftMargin(0.15);
      c1->cd(ip+1)->SetBottomMargin(0.15);
      c1->cd(ip+1)->SetGrid(0,0);
      c1->cd(ip+1)->SetTicks(1,1);
    }
    
    for(int ip = 0; ip < 9; ip++)
    {
      c1->cd(ip+1);
      hist2D_before[ip]->Draw("colz");
    }

    c1->SaveAs("figures/QAbefore.pdf");
 
    TCanvas *c2 = new TCanvas("c2","c2",10,10,900,900);
    c2->Divide(3,3);
    for(int ip = 0; ip < 9; ip++)
    {
      c2->cd(ip+1)->SetLeftMargin(0.15);
      c2->cd(ip+1)->SetBottomMargin(0.15);
      c2->cd(ip+1)->SetGrid(0,0);
      c2->cd(ip+1)->SetTicks(1,1);
    }
    
    for(int ip = 0; ip < 9; ip++)
    {
      c2->cd(ip+1);
      hist2D[ip]->Draw("colz");
    }

    c2->SaveAs("figures/QAafter.pdf");
    std::cout << std::setprecision(10) << lengthMidPin << std::endl;

    TCanvas *c3 = new TCanvas("c3","c3",10,10,300,300);
    c3->cd();
    gRotVLength->SetMarkerStyle(21);
    gRotVLength->Draw("ap");

    c3->SaveAs("figures/QARotPin.pdf");
   
    TCanvas *c4 = new TCanvas("c4","c4",10,10,300,300);
    c4->cd();
    gOriginVLength->SetMarkerStyle(21);
    gOriginVLength->Draw("ap");

    c4->SaveAs("figures/QAOriginPin.pdf");


}
