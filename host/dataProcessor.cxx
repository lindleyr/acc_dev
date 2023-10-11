#include "plotHelper.h"
#include "HoughHelper.h"
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <array>
#include <cmath>
#include <unordered_set>
#include <memory>
#include <chrono>
using namespace std;
using namespace std::chrono;
#define NEVENTS 9996

template <class T> class SmartPtr{

  T* ptr;

public:

  explicit SmartPtr(T* p = NULL) { ptr = p; }

  ~SmartPtr() { delete (ptr); }

  T& operator*() { return *ptr;}

  T* operator->() { return ptr;}

};

int main(int argc,char *argv[]){

  std::string inDir, outDir, file;
  static struct option long_options[] =
  {
    {"inDir", 1, NULL, 'a'},
    {"outDir", 1, NULL, 'b'},
    {"data", 1, NULL, 'c'},
    {NULL, 0, NULL, 0}
  };

  int opt;
  while ( (opt = getopt_long(argc, argv,"abc", long_options, NULL)) != -1 ) {  // for each option...
    switch ( opt )
      {
      case 'a': inDir = optarg; break;
      case 'b': outDir = optarg; break;
      case 'c': file = optarg; break;
      case 0: break;
      }
  }
  // std::vector<float> datavec;
  std::vector<std::vector<float>> datavec(NEVENTS);
  auto start = high_resolution_clock::now();
  GetInfoFromFile_2(file, datavec);
  print_info_vec_data(datavec, 10); // you could do this just to check
  int nlines = 8884; // nlines is different than nevents, each event can have 8+ hits
  
  for(int event=0;event<1;event++){
    int eventnumber = 2;
    unsigned int size_arr = 0;
    std::cout << "size of event " << event << " is " << datavec[event].size() << endl;
    size_arr += datavec[event].size();
    size_arr += 4*datavec[event].size()/9;
    std::cout << "size of array is " << size_arr << endl;
    auto data_arr = std::unique_ptr<double[]>(new double[size_arr]);
    ConvertVecToArr_3(datavec, data_arr, event);
    std::cout << "conversion successful" << endl;
    for(int i=0; i<sizeof(data_arr)/sizeof(double); i++){
       std::cout << "Entry " << i << " of array is " << data_arr[i] << endl;
    }
    //print_info_array_data(data_arr,size_arr);

    //SelectEvent(data_arr, nlines);
    //SelectEvent(data_arr, 1);
    HoughTransformAvg(data_arr,size_arr,event);
    //delete(data_arr);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop-start);
  std::cout << "Time taken is " << duration.count() << " microseconds." << endl;

  return 0;

}
