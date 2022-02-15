#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <numeric>
#include "TMath.h"

void fillMisalignmentTables(std::string date = "20220214")
{
  std::string sensor[3] = {"Inner Sensor (0)", "Outer Sensor (1)", "Outer Sensor (2)"};  

  std::string v_moduleId[108];
  float v_rot[108], v_dx[108], v_dy[108];;
  std::string line;
  std::string moduleId;
  float rot, dx, dy;
  int count = 0;
 
  std::fstream infile_inner;
  infile_inner.open("./InnerSensors/output/InnerAlignmentValues.txt", ios::in);

  if(infile_inner.is_open())
  {
    while(std::getline(infile_inner, line))
    {
      std::stringstream str(line);
      
      if(count == 0)
      {
        count++;
        std::cout << str.str() << endl;
        continue;
      }

      str >> moduleId; v_moduleId[(count-1)*3] = moduleId;        
      str >> rot;      v_rot[(count-1)*3] = rot;            
      str >> dx;       v_dx[(count-1)*3] = dx;
      str >> dy;       v_dy[(count-1)*3] = dy;

      std::cout << count++ << endl; 
      std::cout << moduleId << " " << rot << " " << dx << " " << dy << std::endl;
    }
  }
  
  std::string line1, line2;
  std::fstream infile_outer;
  infile_outer.open("./OuterSensors/output/OuterAlignmentValues.txt", ios::in);

  count = 0;
  if(infile_outer.is_open())
  {
    while(std::getline(infile_outer, line1))
    {
      std::stringstream str1(line1);
      
      if(count == 0)
      {
        count++;
        std::cout << str1.str() << endl;
        continue;
      }

      str1 >> moduleId; v_moduleId[1+(count-1)*3] = moduleId;        
      str1 >> rot;      v_rot[1+(count-1)*3] = rot;            
      str1 >> dx;       v_dx[1+(count-1)*3] = dx;
      str1 >> dy;       v_dy[1+(count-1)*3] = dy;

      std::cout << count << endl; 
      std::cout << moduleId << " " << rot << " " << dx << " " << dy << std::endl;
     
      std::getline(infile_outer, line2);
     
      std::stringstream str2(line2);
      
      str2 >> moduleId; v_moduleId[2+(count-1)*3] = moduleId;        
      str2 >> rot;      v_rot[2+(count-1)*3] = rot;            
      str2 >> dx;       v_dx[2+(count-1)*3] = dx;
      str2 >> dy;       v_dy[2+(count-1)*3] = dy;

      std::cout << count++ << endl; 
      std::cout << moduleId << " " << rot << " " << dx << " " << dy << std::endl;
    }
  }

  std::fstream outfile;
  outfile.open("output/fstSensorOnWedge."+date+".000001.C", ios::out);

  outfile << "TDataSet *CreateTable() {" << std::endl;
  outfile << "    if (!TClass::GetClass(\"St_Survey\")) return 0;" << std::endl;
  outfile << "Survey_st row;" << std::endl;
  outfile << "St_Survey *tableSet = new St_Survey(\"fstSensorOnWedge\",108);" << std::endl;

  for(int i = 0; i < 108; i++)
  {
    outfile <<                                                      std::endl;
    outfile << "memset(&row,0,tableSet->GetRowSize());"          << std::endl;
    outfile << "    row.Id   = " <<  i + 1                << ";" << std::endl;
    outfile << "    row.r00  = " <<  TMath::Cos(v_rot[i]) << ";" << std::endl;
    outfile << "    row.r01  = " << -TMath::Sin(v_rot[i]) << ";" << std::endl;
    outfile << "    row.r02  = 0.0;"                             << std::endl;
    outfile << "    row.r10  = " <<  TMath::Sin(v_rot[i]) << ";" << std::endl;
    outfile << "    row.r11  = " <<  TMath::Cos(v_rot[i]) << ";" << std::endl;
    outfile << "    row.r12  = 0.0;"                             << std::endl; 
    outfile << "    row.r20  = 0.0;"                             << std::endl; 
    outfile << "    row.r21  = 0.0;"                             << std::endl;    
    outfile << "    row.r22  = 1.0;"                             << std::endl;  
    outfile << "    row.t0   = " <<  v_dx[i]              << ";" << std::endl;
    outfile << "    row.t1   = " <<  v_dy[i]              << ";" << std::endl;
    outfile << "    row.t2   = 0.0;"                             << std::endl;
    outfile << "    memcpy(&row.comment,\"Wedge " << int(i/3) << ", " << v_moduleId[i] << ", " << sensor[i%3] << "\\x00\",1);" << std::endl;
    outfile << "tableSet->AddAt(&row);"                          << std::endl;

    std::cout << "    memcpy(&row.comment,\"Wedge " << int(i/3) << ", " << v_moduleId[i] << ", " << sensor[i%3] << "\\x00\",1);" << std::endl;
  }
  
  outfile.close();

}

