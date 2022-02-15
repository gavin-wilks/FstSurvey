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
  std::fstream infile;
  infile.open("data/MisalignmentTables_Wedges.csv", ios::in);

  std::vector<std::string> v_moduleId;
  std::vector<int> v_wedgeId;
  std::vector<float> v_rot, v_dx, v_dy;
  std::string line;
  std::string moduleId;
  int wedgeId;
  float rot, dx, dy;
  char ch;
  int count = 0;
 
  if(infile.is_open())
  {
    while(std::getline(infile, line))
    {
      std::stringstream str(line);
      
      if(count == 0)
      {
        count++;
        std::cout << str.str() << endl;
        continue;
      }

      std::getline(str, moduleId, ','); v_moduleId.push_back(moduleId); 
      str >> wedgeId;        str >> ch; v_wedgeId.push_back(wedgeId);  
      str >> rot;            str >> ch; v_rot.push_back(rot);           
      str >> dx;             str >> ch; v_dx.push_back(dx/10.0);             
      str >> dy;                        v_dy.push_back(dy/10.0);             

      std::cout << count++ << endl; 
      std::cout << moduleId << "," << wedgeId << "," << rot << "," << dx << "," << dy << std::endl;
    }
  }

  vector<size_t> idx(v_wedgeId.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v_wedgeId](size_t i1, size_t i2) {return v_wedgeId[i1] < v_wedgeId[i2];});

  std::fstream outfile;
  outfile.open("output/fstWedgeOnHss."+date+".000001.C", ios::out);

  outfile << "TDataSet *CreateTable() {" << std::endl;
  outfile << "    if (!TClass::GetClass(\"St_Survey\")) return 0;" << std::endl;
  outfile << "Survey_st row;" << std::endl;
  outfile << "St_Survey *tableSet = new St_Survey(\"fstWedgeOnHss\",36);" << std::endl;

  for( auto i : idx)
  {
    outfile <<                                                      std::endl;
    outfile << "memset(&row,0,tableSet->GetRowSize());"          << std::endl;
    outfile << "    row.Id   = " <<  v_wedgeId[i] + 1     << ";" << std::endl;
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
    outfile << "    memcpy(&row.comment,\"Wedge " << v_wedgeId[i] << ", " << v_moduleId[i] << "\\x00\",1);" << std::endl;
    outfile << "tableSet->AddAt(&row);"                          << std::endl;
  }
  
  outfile.close();

}

