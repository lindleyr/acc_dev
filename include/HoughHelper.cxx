#include <dirent.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include "HoughHelper.h"
#include <memory>
#include <bits/stdc++.h>
#ifndef HoughHelper_cxx
#define HoughHelper_cxx

using namespace std;

// ================================================
// ================================================
// Structs to read the txt files


pvec operator+(const pvec& a, const pvec& b) {
    return { {a[0]+b[0], a[1]+b[1]} };
}
pvec operator-(const pvec& a, const pvec& b) {
    return { {a[0]-b[0], a[1]-b[1]} };
}
pvec operator*(const pvec& a, double scale) {
    return { {a[0]*scale, a[1]*scale} };
}

double length(const pvec& v) {
    return std::hypot(v[0], v[1]);
}

pvec rotate90( const pvec& v) {
    return { -v[1], v[0]};
}

double crossProduct( const pvec& a, const pvec& b ) {
    return a[0]*b[1] - a[1]*b[0];
}


// ================================================
// ================================================
void GetInfoFromFiles(string hitsFile, string partsFile, hit *hits, particle *particles, int nevents){
  // initialize structs (maybe cleaner way?)
  cout << " i have started initialize " << endl;

  for (int i = 0; i < nevents; ++i){
    for (int j = 0; j < 100; ++j){
      hits[i].layer[j] = -9999; hits[i].r[j] = -9999; hits[i].x[j] = -9999; hits[i].y[j] = -9999; hits[i].z[j] = -9999;
      particles[i].barcode[j] = -9999; particles[i].charge[j] = -9999; particles[i].pt[j] = -9999; particles[i].d0[j] = -9999;
    }
    /*for (int j = 0; j < 100; ++j){                                                                                                   
      hits[i].layer[j] = -9999; hits[i].r[j] = -9999; hits[i].x[j] = -9999; hits[i].y[j] = -9999; hits[i].z[j] = -9999;                         particles[i].charge[j] = -9999; particles[i].pt[j] = -9999; particles[i].d0[j] = -9999;
    }*/
  }

  cout << " i have finished initialize " << endl;

  // hits file
  int event;
  double layer, r, x, y, z;
  std::string line;
  std::ifstream hitNameFile(hitsFile.c_str());
  cout << " hits filename : " << hitsFile << endl;
  int event1=0;
  int idx=-1;
  while (std::getline(hitNameFile, line)){
    std::stringstream ss(line);
    ss >> event >> layer >> r >> x >> y >> z;
    if( event1 == event){ idx++; } else { idx = 0; }
    hits[event].layer[idx] = layer; hits[event].r[idx] = r; hits[event].x[idx] = x; hits[event].y[idx] = y; hits[event].z[idx] = z;
    event1 = event;
  }

  // particles file
  int event_part;
  double barcode, charge, pt, d0;
  std::string line_part;
  std::ifstream partNameFile(partsFile.c_str());
  cout << " particles filename : " << partsFile << endl;
  int event1_part=0;
  int idx_part=-1;
  while (std::getline(partNameFile, line_part)){
    std::stringstream ss(line_part);
    ss >> event_part >> barcode >> charge >> pt >> d0 ;
    if( event1_part == event_part){ idx_part++; } else { idx_part = 0; }
      particles[event_part].barcode[idx_part] = barcode; particles[event_part].charge[idx_part] = charge; particles[event_part].pt[idx_part] = pt; particles[event_part].d0[idx_part] = d0;
  }
}

void GetInfoFromFile(string mergeFile, std::vector<float>& vec){
  cout << " i am in: GetInfoFromFile " << endl;

  int event;
  double layer, r, x, y, z, charge, pt, d0, numhits;
  std::string line;
  std::ifstream MergeNameFile(mergeFile.c_str());
  cout << " data filename : " << mergeFile << endl;
  while (std::getline(MergeNameFile, line)){
    std::stringstream ss(line);
    ss >> event >> layer >> r >> x >> y >> z >> charge >> pt >> d0 >> numhits;
    vec.push_back(event);
    vec.push_back(layer);
    vec.push_back(r);
    vec.push_back(x);
    vec.push_back(y);
    vec.push_back(z);
    vec.push_back(charge);
    vec.push_back(pt);
    vec.push_back(d0);
    vec.push_back(numhits);
  }
}

void GetInfoFromFile_2(string mergeFile, std::vector<std::vector<float>>& vec){
  cout << " i am in: GetInfoFromFile " << endl;

  int event;
  double layer, r, x, y, z, barcode, charge, pt, d0;
  long double hbarcode;
  std::string line;
  cout << "line max size " << line.max_size() << endl;
  std::ifstream MergeNameFile(mergeFile.c_str());
  cout << " data filename : " << mergeFile << endl;
  // int cache=-1;                                                                                                                                                                                     
  while (std::getline(MergeNameFile, line)){
    std::stringstream ss(line);
    ss >> event >> layer >> r >> x >> y >> z >> hbarcode >> barcode >> charge >> pt >> d0;
    //cout << event << " | " << layer << " | " << r << " | " << x << " | " << y << " | " << " | " << z << " | " << hbarcode << " | " << barcode << " | " << charge << " | " << pt << " | " << d0 << endl;

    vec[event].push_back(layer); //cout << "got layer" << endl;
    vec[event].push_back(r); //cout << "got r" << endl;
    vec[event].push_back(x); //cout << "got x" << endl;
    vec[event].push_back(y);
    vec[event].push_back(z);
    vec[event].push_back(hbarcode);
    vec[event].push_back(barcode);
    vec[event].push_back(charge);
    vec[event].push_back(pt);
    vec[event].push_back(d0);
    // cache=event;                                                                                                                                                                                    
  }
}


void GetInfoFromFile(string mergeFile, std::vector<std::vector<float>>& vec){
  cout << " i am in: GetInfoFromFile " << endl;

  int event;
  double layer, r, x, y, z, charge, pt, d0, numhits;
  std::string line;
  std::ifstream MergeNameFile(mergeFile.c_str());
  cout << " data filename : " << mergeFile << endl;
  // int cache=-1;
  while (std::getline(MergeNameFile, line)){
    std::stringstream ss(line);
    ss >> event >> layer >> r >> x >> y >> z >> charge >> pt >> d0 >> numhits;
    // if (event>cache) vec[event].push_back(numhits);
    vec[event].push_back(numhits);
    vec[event].push_back(layer);
    vec[event].push_back(r);
    vec[event].push_back(x);
    vec[event].push_back(y);
    vec[event].push_back(z);
    vec[event].push_back(charge);
    vec[event].push_back(pt);
    vec[event].push_back(d0);
    // cache=event;
  }
}

void print_info_vec_data(std::vector<float>& vec, int size){

  for(int i=0; i<size; ++i){
    std::cout << " event "<< vec[10*i]
              << " layer " << vec[10*i+1]
              << " r " << vec[10*i+2]
              << " x " << vec[10*i+3]
              << " y " << vec[10*i+4]
              << " z " << vec[10*i+5]
              << " charge " << vec[10*i+6]
              << " pt " << vec[10*i+7]
              << " d0 " << vec[10*i+8]
              << " numhits " << vec[10*i+9]
              << std::endl;
  }
}

void print_info_vec_data(std::vector<std::vector<float>>& vec, int size){

  for(unsigned int i=0; i<size; i++){
    std::cout << " eventNumber: #" << i<< std::endl;
    for(unsigned int j=0; j<vec[i].size(); j++){
      std::cout << " " <<vec[i][j]<< " ";
    }
    std::cout << std::endl;
    std::cout << " ============================== " << std::endl;
  }
}

//void ConvertVecToArr_3(std::vector<std::vector<float>>& vec, std::unique_ptr<double[]> *arr, int eventnumber){

void ConvertVecToArr_3(std::vector<std::vector<float>>& vec, std::unique_ptr<double[]>& arr, int eventnumber){
  int numhits = 0;
  if(vec[eventnumber].size() > 0){
    numhits = vec[eventnumber].size()/10;
    std::cout << "eventnumber is " << eventnumber << " and numhits is " << numhits << endl;
    for(int j=0; j<numhits; j++){
      arr[14*(j)]=eventnumber; //eventnumber
      arr[14*(j)+1]=j;         // hit number
      arr[14*(j)+2]=vec[eventnumber][j*10]; //layer
      arr[14*(j)+3]=vec[eventnumber][j*10+1];  //r
      arr[14*(j)+4]=vec[eventnumber][j*10+2];  //x
      arr[14*(j)+5]=vec[eventnumber][j*10+3];  //y
      arr[14*(j)+6]=vec[eventnumber][j*10+4];  //z
      arr[14*(j)+7]=vec[eventnumber][j*10+5];  //hbarcode
      arr[14*(j)+8]=vec[eventnumber][j*10+6];  //barcode
      arr[14*(j)+9]=vec[eventnumber][j*10+7];  //charge
      arr[14*(j)+10]=vec[eventnumber][j*10+8];  //pT
      arr[14*(j)+11]=vec[eventnumber][j*10+9]; //d0 
      arr[14*(j)+12]=0; // # roads contributed to
      arr[14*(j)+13]=0; // total # roads
    }

  }

}

void ConvertVecToArr_2(std::vector<std::vector<float>>& vec, double *arr, int size){
  int numhits = 0;
  for (unsigned int i = 0; i < size; i++) {
    if(vec[i].size() > 0){
      numhits = vec[i].size()/9;
      std::cout << "i is " << i << " and numhits is " << numhits << endl;                                                                                                                           
      for(int j=0; j<numhits; j++){
	//std::cout << "j is " << j << endl;                                                                                                                                                          
	arr[10*(i+j)]=i;
	//std::cout << "event is " << arr[10*(i+j)] << " and should be " << i;                                                                                                                        
	//arr[10*(i+j)+1]=vec[i][j*9];
	//std::cout << " numhits is " << arr[10*(i+j)+1] << " and should be " << vec[i][j*9] ;                                                                                                        
	arr[10*(i+j)+1]=vec[i][j*9];
	//std::cout << " layer " << arr[10*(i+j)+1];                                                                                                                                                  
	arr[10*(i+j)+2]=vec[i][j*9+1];
	//std::cout << " r " << arr[10*(i+j)+2];
	arr[10*(i+j)+3]=vec[i][j*9+2];
	//std::cout << " x " << arr[10*(i+j)+3];                                                                                                                                                      
	arr[10*(i+j)+4]=vec[i][j*9+3];
	//std::cout << " y " << arr[10*(i+j)+4];                                                                                                                                                      
	arr[10*(i+j)+5]=vec[i][j*9+4];
	//std::cout << " z " << arr[10*(i+j)+5];                                                                                                                                                      
	arr[10*(i+j)+6]=vec[i][j*9+5];
	//std::cout << "barcode " << arr[10*(i+j)+6];
        arr[10*(i+j)+7]=vec[i][j*9+6];
	//std::cout << " charge " << arr[10*(i+j)+7];                                                                                                                                                 
	arr[10*(i+j)+8]=vec[i][j*9+7];
	//std::cout << " pt " << arr[10*(i+j)+8];                                                                                                                                                     
	arr[10*(i+j)+9]=vec[i][j*9+8];
	//std::cout << " d0 " << arr[10*(i+j)+9];                                                                                                                                                     
	//std::cout << endl;                                                                                                                                                                          
      }
    }
  }

}


void ConvertVecToArr(std::vector<std::vector<float>>& vec, double *arr, int size){
     int numhits = 0;
     for (unsigned int i = 0; i < size; i++) {
       if(vec[i].size() > 0){
       numhits = vec[i][0];
       //std::cout << "i is " << i << " and numhits is " << numhits << endl;
       for(int j=0; j<numhits; j++){
	 //std::cout << "j is " << j << endl;
         arr[10*(i+j)]=i;
         //std::cout << "event is " << arr[10*(i+j)] << " and should be " << i;
         arr[10*(i+j)+1]=vec[i][j*9];
         //std::cout << " numhits is " << arr[10*(i+j)+1] << " and should be " << vec[i][j*9] ;
         arr[10*(i+j)+2]=vec[i][j*9+1];
         //std::cout << " layer " << arr[10*(i+j)+2];
         arr[10*(i+j)+3]=vec[i][j*9+2];
         //std::cout << " r " << arr[10*(i+j)+3];
         arr[10*(i+j)+4]=vec[i][j*9+3];
         //std::cout << " x " << arr[10*(i+j)+4];
         arr[10*(i+j)+5]=vec[i][j*9+4];
         //std::cout << " y " << arr[10*(i+j)+5];
         arr[10*(i+j)+6]=vec[i][j*9+5];
         //std::cout << " z " << arr[10*(i+j)+6];
         arr[10*(i+j)+7]=vec[i][j*9+6];
         //std::cout << " charge " << arr[10*(i+j)+7];
         arr[10*(i+j)+8]=vec[i][j*9+7];
         //std::cout << " pt " << arr[10*(i+j)+8];
         arr[10*(i+j)+9]=vec[i][j*9+8];
         //std::cout << " d0 " << arr[10*(i+j)+9];
         //std::cout << endl;
       }
       }
     } 

}

//
// void print_info_array_data(double *arr, unsigned int size){
//
//   for (unsigned int i = 0; i < size; i++) {
//     std::cout << " event "<< arr[10*i]
//               << " layer " << arr[10*i+1]
//               << " r " << arr[10*i+2]
//               << " x " << arr[10*i+3]
//               << " y " << arr[10*i+4]
//               << " z " << arr[10*i+5]
//               << " charge " << arr[10*i+6]
//               << " pt " << arr[10*i+7]
//               << " d0 " << arr[10*i+8]
//               << " numhits " << arr[10*i+9]
//               << std::endl;
//   }
// }

//void print_info_array_data(double *arr, unsigned int size){
void print_info_array_data(std::unique_ptr<double[]>& arr, int arrsize){ //focus on one event for now
  // numhits layer r x y z charge pt d0
  //int arrsize = *(&arr+1)-arr;
  std::cout << "size of array is " << arrsize << endl;
  for (unsigned int i = 0; i < arrsize/14; i++) {
    std::cout << " event " << arr[14*i]
              << " hit number "<< arr[14*i+1]
              << " layer " << arr[14*i+2]
              << " r " << arr[14*i+3]
              << " x " << arr[14*i+4]
              << " y " << arr[14*i+5]
              << " z " << arr[14*i+6]
              << " hbarcode " << arr[14*i+7]
              << " barcode " << arr[14*i+8]
              << " charge " << arr[14*i+9]
              << " pt " << arr[14*i+10]
              << " d0 " << arr[14*i+11]
              << std::endl;
  }
}


void SelectEvent(double *arr, unsigned int size){
  //
  // const double m_acceptedDistanceBetweenLayersMin = 200; // min R disstance for hits pair filtering
  // const double m_acceptedDistanceBetweenLayersMax = 600;
  //
  // float m_d0_range = 120;
  // float m_qOverPt_range = 0.002;
  // int m_imageSize_x = 216; // i.e. number of bins in d0
  // int m_imageSize_y = 216; // i.e. number of bins in q/pT
  // double m_step_x = (2*m_d0_range) / m_imageSize_x; // helpers (accumulator granularity)
  // double m_step_y = (2*m_qOverPt_range) / m_imageSize_y;
  // bool m_continuous = true; // assure that there is continuity of the line (i.e. middle bins in d0 are filled when one q/pT step would result in a hole)

  // Loop over size of array
  for (unsigned int i = 0; i <size; i++) {
    std::cout << " i = "<< i
              << " event "<< arr[10*i]
              << " layer " << arr[10*i+1]
              << " r " << arr[10*i+2]
              << " x " << arr[10*i+3]
              << " y " << arr[10*i+4]
              << " z " << arr[10*i+5]
              << " charge " << arr[10*i+6]
              << " pt " << arr[10*i+7]
              << " d0 " << arr[10*i+8]
              << " numhits " << arr[10*i+9]
              << std::endl;
  }
  // check if events have the same eventNumber
  // Loop over hit 1
  // Loop over hit 2
  // Put conditions

}



// ================================================
// ===================================================
void printout_struct(hit *hits, particle *particles){
  // Check outputs
  cout << particles[1].barcode[0] << " " << particles[1].charge[0] << " " << particles[1].pt[0]  << " " << particles[1].d0[0] << endl;
  for(int ii=0; ii<100; ii++){
    if (hits[0].layer[ii]!=-9999){
      cout << hits[0].layer[ii] << " " << hits[0].r[ii] << " " << hits[0].x[ii]  << " " << hits[0].y[ii] << " " << hits[0].z[ii] << endl;
    } else {
      continue;
    }
  }
}
// ================================================
// ===================================================
// bool passThreshold(vector2D<std::pair<int, vector<double>>> &image, int x, int y,  double m_step_x, double m_d0_range) {
bool passThreshold(vector2D<std::pair<int, hit>> &image, int x, int y,  double m_step_x, double m_d0_range, hit *hits, particle *particles) {
    const int count = image(x,y).first;
    const float d0 = xtod0(x, m_step_x, m_d0_range);
    std::cout << "count is " << count << " and d0 is " << d0 << endl;
    int m_threshold = 8;
    int m_threshold50 = 8;
    if ( std::abs(d0) < 50.0 && count >= m_threshold50 ){
      std::cout << "Event passes threshold for d0 below 50!" << endl;
      return true;
    }
    if ( std::abs(d0) >= 50.0 && count >= m_threshold ){
      std::cout << "Event passes threshold for d0 above 50!" << endl;
      return true;
    }
    return false;
}
bool passThreshold(vector2D<std::pair<int, hit>> &image, int x, int y,  double m_step_x, double m_d0_range, int numhits) {
    const int count = image(x,y).first;
    const float d0 = xtod0(x, m_step_x, m_d0_range);
    //std::cout << "d0 for (89.5,98.5) is " << xtod0(89.5, m_step_x, m_d0_range) << endl; 
    int m_threshold = 15;
    int m_threshold50 = 15;
    if ( std::abs(d0) < 50.0 && count >= m_threshold50 ){
      //std::cout << "d0 is " << d0 << " and count is " << count << ". Event passes threshold for d0 below 50!" << endl;
      return true;
    }
    if ( std::abs(d0) >= 50.0 && count >= m_threshold ){
      //std::cout << "d0 is " << d0 << " and count is " << count << ". Event passes threshold for d0 below 50!" << endl;
      return true;
    }
    return false;
}

// bool isLocalMaxima(vector2D<std::pair<int, vector<double>>> &image, int x, int y, int m_imageSize_x, int m_imageSize_y) {
bool isLocalMaxima(vector2D<std::pair<int, hit>> &image, int x, int y, int m_imageSize_x, int m_imageSize_y, hit *hits, particle *particles) {
    const auto centerValue =  image(x,y).first;
    for ( int xaround = std::min(x-1, 0); xaround <= std::max(m_imageSize_x-1, x+1); xaround++  ) {
        for ( int yaround = std::min(y-1, 0); yaround <= std::max(m_imageSize_y-1, y+1); yaround++  ) {
        if ( image(xaround,yaround).first > centerValue ) { return false; }
        }
    }
    return true;
}
bool isLocalMaxima(vector2D<std::pair<int, hit>> &image, int x, int y, int m_imageSize_x, int m_imageSize_y) {
   const auto centerValue =  image(x,y).first;
   for ( int xaround = std::max(x-2, 0); xaround <= std::min(m_imageSize_x-2, x+2); xaround++  ) {
     for ( int yaround = std::max(y-2, 0); yaround <= std::min(m_imageSize_y-2, y+2); yaround++  ) {
       //std::cout << "xaround is " << xaround << " and yaround is " << yaround << endl; 
       // for ( int xaround = std::min(x-1, 0); xaround <= std::max(m_imageSize_x-1, x+1); xaround++  ) {
       // for ( int yaround = std::min(y-1, 0); yaround <= std::max(m_imageSize_y-1, y+1); yaround++  ) {
       //       if(xaround != x || yaround != y){
       //std::cout << "centerValue for x " << x << " and y " << y << " is " << centerValue << " and for xaround " << xaround << " and yaround " << yaround << " is " << image(xaround,yaround).first << endl;  
         if ( image(xaround,yaround).first > centerValue) { return false; }
	 //}
     }  
   }
    return true;
}

pvec LocalMaximaTest(vector2D<std::pair<int, hit>> &image, int x, int y, int m_imageSize_x, int m_imageSize_y, double m_step_x, double m_step_y, double m_d0_range, double m_qOverPt_range) {
  const auto centerValue =  image(x,y).first;
  std::cout << "x is " << x << ", y is " << y << ", count is " << centerValue << endl;
  for ( int xaround = std::max(x-2, 0); xaround <= std::min(m_imageSize_x-2, x+2); xaround++  ) {
    for ( int yaround = std::max(y-2, 0); yaround <= std::min(m_imageSize_y-2, y+2); yaround++  ) {
      //std::cout << "xaround is " << xaround << " and yaround is " << yaround << endl;
      if ( image(xaround,yaround).first > centerValue) { return {{-1,-1}}; }
      if ( image(xaround,yaround).first == centerValue && (xaround != x || yaround != y)){
	const float d0avg = 0.5 * (xtod0(x, m_step_x, m_d0_range) + xtod0(xaround, m_step_x, m_d0_range));
	const float qptavg = 0.5 * (ytoqoverpt(y, m_step_y, m_qOverPt_range) + ytoqoverpt(yaround, m_step_y, m_qOverPt_range));
	if((yaround*215+xaround) > (y*215+x)){
	  return {{-1, -1}};
        }
        if((yaround*215+xaround) < (y*215+x)){
	  return {{d0avg, qptavg}};
	}
      }
    }
  }
  const float d0 = xtod0(x, m_step_x, m_d0_range);
  //std::cout << "d0 is " << d0 << endl;
  return {{d0,ytoqoverpt(y, m_step_y, m_qOverPt_range)}};
}



void SelectEvents(hit *hits, particle *particles, int nevents){

  const double m_acceptedDistanceBetweenLayersMin = 200; // min R disstance for hits pair filtering
  const double m_acceptedDistanceBetweenLayersMax = 600;

  float m_d0_range = 120;
  float m_qOverPt_range = 0.002;
  int m_imageSize_x = 216; // i.e. number of bins in d0
  int m_imageSize_y = 216; // i.e. number of bins in q/pT
  double m_step_x = (2*m_d0_range) / m_imageSize_x; // helpers (accumulator granularity)
  double m_step_y = (2*m_qOverPt_range) / m_imageSize_y;
  bool m_continuous = true; // assure that there is continuity of the line (i.e. middle bins in d0 are filled when one q/pT step would result in a hole)

  for(int event=0; event<nevents; event++){
    // if (event>2) continue;
    vector2D<std::pair<int, hit>> image(m_imageSize_x, m_imageSize_y);

    for(int ihit1=0; ihit1<100; ihit1++){
      for(int ihit2=ihit1+1; ihit2<100; ihit2++){

        double radius_hit1 = GetR(hits[event].x[ihit1], hits[event].y[ihit1]);
        double radius_hit2 = GetR(hits[event].x[ihit2], hits[event].y[ihit2]);
        double radiusDifference =  radius_hit2 - radius_hit1;

        if ( hits[event].layer[ihit1] == hits[event].layer[ihit2]){
          continue;
        }

        if (  not (m_acceptedDistanceBetweenLayersMin < radiusDifference && radiusDifference < m_acceptedDistanceBetweenLayersMax) ){
          continue;
        }

        const pvec p1 {{hits[event].x[ihit1], hits[event].y[ihit1]}};
        const pvec p2 {{hits[event].x[ihit2], hits[event].y[ihit2]}};
        const pvec halfDiff = (p2 - p1)*0.5;
        const fp_t halfLen = length(halfDiff);

        int xbefore = -1;

        for ( int y = 1; y < m_imageSize_y+1; y++ ) {
          const fp_t qoverpt = -1.*( (y * m_step_y) + m_step_y*0.5 - m_qOverPt_range);
          const fp_t radius = 1.0/(0.6*qoverpt);
          const fp_t scale = std::copysign( std::sqrt( std::pow(radius/halfLen, 2) - 1), radius );
          const pvec rprime = rotate90(halfDiff) * scale;
          const pvec center = p1 + halfDiff + rprime;
          const fp_t d0 =  (std::signbit(radius) ? -1.0 : 1.0)*(length(center) - abs(radius));
          int x = (d0 + m_d0_range) / m_step_x;
          // cout << x << endl;
          if ( 1 <= x && x < m_imageSize_x) {
            if (xbefore == -1) xbefore = x;
            if ( m_continuous ) { // fill the bins along x starting from the last one filled
              const int xmin =  (xbefore < x)? xbefore: x;
              const int xmax =  (xbefore < x)? x: xbefore;
              for ( int xinterpolated = xmin; xinterpolated <= xmax; ++xinterpolated) {
                image(xinterpolated, y).first++;
                image(xinterpolated, y).second =  hits[event] ;
                // image(xinterpolated, y).second =  hits[event].x[ihit2] ;
              }
            } else {
                image(x, y).first++;
                image(x, y).second =  hits[event] ;
                // image(x, y).second.insert( hits[event].x );
                // image(x, y).second.insert( hits[event].x );
            }
            xbefore = x;
            }
          }
      }
    }


    for (int y = 0; y < m_imageSize_y; y++) {
      for (int x = 0; x < m_imageSize_x; x++) {
        if (passThreshold(image, x, y, m_step_x, m_d0_range, hits, particles) && isLocalMaxima( image, x, y, m_imageSize_x, m_imageSize_y, hits, particles) ) {
          // roads.push_back(createRoad(image(x, y).second, x, y));
          if (event==1 || event==2 || event == 3) {
            // cout << image(x, y).second.x << " " << image(x, y).second.r ;
            cout << " d0: " << xtod0(x, m_step_x, m_d0_range) << " truthd0: " << particles[event].d0[0]  << " resolution d0 :" << (particles[event].d0[0] - xtod0(x, m_step_x, m_d0_range) )
                 << " q/pt " << ytoqoverpt(y, m_step_y, m_qOverPt_range) << " truth q/pT: " << particles[event].charge[0] / particles[event].pt[0]<< " resolution q/pT :" << (particles[event].charge[0] / particles[event].pt[0] - ytoqoverpt(y, m_step_y, m_qOverPt_range))<<  endl;
          }
        }
      }
    }
  }
}


void HoughTransform(std::unique_ptr<double[]>& arr, int arrsize, int eventnumber){


  std::cout << " Hello I am here " << std::endl;
  float m_d0_range = 120;
  float m_qOverPt_range = 0.002;
  int m_imageSize_x = 216; // i.e. number of bins in d0 (nominal 216)
  int m_imageSize_y = 216; // i.e. number of bins in q/pT (nominal 216)
  double m_step_x = (2*m_d0_range) / m_imageSize_x; // helpers (accumulator granularity)
  double m_step_y = (2*m_qOverPt_range) / m_imageSize_y;
  //std::cout << "m_step_x is " << m_step_x << " and m_step_y is " << m_step_y << endl;
  bool m_continuous = true; // assure that there is continuity of the line (i.e. middle bins in d0 are filled when one q/pT step would result in a hole)

  // if (event>2) continue;
  vector2D<std::pair<int, hit>> image(m_imageSize_x, m_imageSize_y);

  //int sizeout = sizeof(arr);
  int sizeout = arrsize/14;
  // sizeout /= 9;
  std::cout << " sizeout: " << sizeout<<  std::endl;

  std::string plotname = "outputfiles/";
  plotname.append(std::to_string(eventnumber));
  plotname.append("_outputvec_nopileup_withOuterLayers.txt");
  char nplot[plotname.size()+1];
  strcpy(nplot,plotname.c_str());
  ofstream myfile;
  myfile.open(nplot);
  //std::cout << nplot << " has been opened." << endl;
  std::vector<std::vector<int>> pix_list(m_imageSize_x*m_imageSize_y);
  int pixel = 0;

  for(int ihit1=0; ihit1<sizeout; ihit1++){
    for(int ihit2=ihit1+1; ihit2<sizeout; ihit2++){

      const pvec p1 {{arr[14*ihit1+4], arr[14*ihit1+5]}}; // x and y for hit1
      //std::cout << "x for hit 1 of " << arr[14*ihit1+1] << " is " << arr[14*ihit1+4] << " and y is " << arr[14*ihit1+5] << endl;
      const pvec p2 {{arr[14*ihit2+4], arr[14*ihit2+5]}};// x and y for hit1
      //std::cout << "x for hit 2 of " << arr[14*ihit2+1] << " is " << arr[14*ihit2+4] << " and y is " << arr[14*ihit2+5] << endl;
      const pvec halfDiff = (p2 - p1)*0.5;
      const fp_t halfLen = length(halfDiff);

      int xbefore = -1;

      for ( int y = 1; y < m_imageSize_y+1; y++ ) {
        const fp_t qoverpt = -1.*( (y * m_step_y) + m_step_y*0.5 - m_qOverPt_range);
	//std::cout << "for y " << y << " qoverpt is " << qoverpt << endl; 
        const fp_t radius = 1.0/(0.6*qoverpt);
	//std::cout << "radius is " << radius << endl;
        const fp_t scale = std::copysign( std::sqrt( std::pow(radius/halfLen, 2) - 1), radius );
	//std::cout << "scale is " << scale << endl;
        const pvec rprime = rotate90(halfDiff) * scale;
        const pvec center = p1 + halfDiff + rprime;
	//std::cout << "p1[0] is " << p1[0] << " and p1[1] is " << p1[1] << endl;
	//std::cout << "halfDiff[0] is " << halfDiff[0] << " and halfDiff[1] is " << halfDiff[1] << endl;
	//std::cout << "rprime[0] is " << rprime[0] << " and rprime[1] is " << rprime[1] << endl;
	//std::cout << "center[0] is " << center[0] << " and center[1] is " << center[1] << endl;
	//std::cout << "length(center) is " << length(center) << endl;
        const fp_t d0 =  (std::signbit(radius) ? -1.0 : 1.0)*(length(center) - abs(radius));
	//		std::cout << "d0 is " << d0 << endl; 
        int x = (d0 + m_d0_range) / m_step_x;
	//        cout << " x: " << x << endl;
        if ( 1 <= x && x < m_imageSize_x) {
          if (xbefore == -1) xbefore = x;
          if ( m_continuous ) { // fill the bins along x starting from the last one filled
            const int xmin =  (xbefore < x)? xbefore: x;
	    //std::cout << "xmin is " << xmin << endl;
            const int xmax =  (xbefore < x)? x: xbefore;
            for ( int xinterpolated = xmin; xinterpolated <= xmax; ++xinterpolated) {
              pixel = (y-1)*m_imageSize_x + xinterpolated;
	      //std::cout << "y is " << y << " x interpolated is " << xinterpolated << " pixel is " << pixel << " hit1 is " << arr[13*ihit1+1] << " hit2 is " << arr[13*ihit2+1] << " max pixlist size is " << m_imageSize_x*m_imageSize_y << endl;
              pix_list[pixel].push_back(arr[14*ihit1+1]);
              pix_list[pixel].push_back(arr[14*ihit2+1]);
	      //std::cout << "added to pixlist" << endl;
              image(xinterpolated, y).first++;
	      myfile << xinterpolated << "," << y-1 << "," << image(xinterpolated,y).first << "," << arr[14*ihit1+1] << "," << arr[14*ihit2+1];
	      myfile << endl;
              // image(xinterpolated, y).second =  hits[event] ;
              // image(xinterpolated, y).second =  hits[event].x[ihit2] ;
              pixel=0;
            }
          } else {
	    //std::cout << "x is " << x << endl;
	      pixel = (y-1)*m_imageSize_x + x;
	      pix_list[pixel].push_back(ihit1);
              pix_list[pixel].push_back(ihit2);
              image(x, y).first++;
	      myfile << x << "," << y-1 << "," << image(x,y).first << "," << arr[14*ihit1+1] << "," << arr[14*ihit2+1];
	      myfile << endl;
              // image(x, y).second =  hits[event] ;
              // image(x, y).second.insert( hits[event].x );
              // image(x, y).second.insert( hits[event].x );
          }
          xbefore = x;
        }
      }
    }
  }

  //sort into four categories:
  //1. Events with 1 hit per layer (Type 0)
  //2. Events with missing layers only (Type 1)
  //3. Events with extra hits only (Type 2)
  //4. Events with both missing layers and extra hits (Type 3)

  int Type = 0;
  double layerarr[sizeout];
  for(int hitties=0; hitties<sizeout; hitties++){
    layerarr[hitties] = arr[14*hitties+2];
  }
  std::sort(layerarr,layerarr+sizeout);
  if(layerarr[0] != 0){Type = 1;}
  for(int i=1; i<sizeout; i++){
    if(layerarr[i-1]+1 == layerarr[i]){
      continue;
    }
    if((layerarr[i-1]+1 > layerarr[i]) && Type == 0){
      Type = 2;
    }
    if((layerarr[i-1]+1 > layerarr[i]) && Type == 1){
      Type = 3;
    }
    if((layerarr[i-1]+1 < layerarr[i]) && Type == 0){
      Type = 1;
    }
    if((layerarr[i-1]+1 < layerarr[i]) && Type == 2){
      Type = 3;
    }
    if(Type == 3) break;
  }
  //std::cout << "Type is " << Type << endl;  
  int nroads = 0;
  int nb50 = 0;
  for (int y = 0; y < m_imageSize_y; y++) {
    for (int x = 0; x < m_imageSize_x; x++) {
      float nmatchedhits = 0;
      int pix = y*m_imageSize_x + x;
      //pvec testresult = LocalMaximaTest(image,x,y,m_imageSize_x,m_imageSize_y,m_step_x,m_step_y,m_d0_range,m_qOverPt_range);
      if (passThreshold(image, x, y, m_step_x, m_d0_range, sizeout) && isLocalMaxima( image, x, y, m_imageSize_x, m_imageSize_y) ) {
      //if (passThreshold(image, x, y, m_step_x, m_d0_range, sizeout)){
	pvec testresult = LocalMaximaTest(image,x,y,m_imageSize_x,m_imageSize_y,m_step_x,m_step_y,m_d0_range,m_qOverPt_range);
	//std::cout << "testresult[x] is " << testresult[0] << endl;
	//if(testresult[0] != -1){
	  nroads++;
	  std::cout << "testresult[x] is " << testresult[0] << " and testresult[y] is " << testresult[1] << endl;
	  cout << " d0: " << xtod0(x, m_step_x, m_d0_range) << " truthd0: " << arr[11]  << " resolution d0 :" << (arr[11] - xtod0(x, m_step_x, m_d0_range) ) << " q/pt " << ytoqoverpt(y,m_step_y,m_qOverPt_range) << " truth q/pT: " << arr[9] / arr[10]<< " resolution q/pT :" << (arr[9] / arr[10] - ytoqoverpt(y,m_step_y,m_qOverPt_range)) <<  endl;
	  for(int i=0; i<pix_list[pix].size(); i++){
	    std::cout << pix_list[pix][i] << ",";
	  }
	
	  //std::cout << endl;
	  for(int j=0; j<sizeout; j++){
	    for(int i=0; i<pix_list[pix].size(); i++){
	      if(j == pix_list[pix][i]){
		//arr[14*j+12]++;
		if(arr[14*j+7] == arr[14*j+8]){ //if hbarcode matches track barcode
		  nmatchedhits++;
		}
	      }
	    }
	  }
	  std::cout << "Barcode fraction for road is " << nmatchedhits/pix_list[pix].size() << endl;
	  if(nmatchedhits/pix_list[pix].size() > 0.5){
	    //arr[14*i+12] = nmatchedhits/pix_list[pix].size();
	    nb50++;
	  } 
	  //}
      }
    }
  }

  ofstream myfile2;
  myfile2.open("outputfiles/outputarr_nopileup_withOuterLayers.txt",std::ios_base::app);

  for(int i=0; i<sizeout; i++){
    arr[14*i+13]=nroads;
    //std::cout << "Efficiency for hit " << i << " is " << arr[14*i+12]/arr[14*i+13] << endl;
    myfile2 << arr[14*i] << " " << arr[14*i+1] << " " << arr[14*i+2] << " " << arr[14*i+3] << " " << arr[14*i+4] << " " << arr[14*i+5] << " " << arr[14*i+6] << " " << arr[14*i+7] << " " << arr[14*i+8] << " " << arr[14*i+9] << " " << arr[14*i+10] << " " << arr[14*i+11] << " " << nb50 << " " << arr[14*i+13] << " " << Type << endl;				       }

}

void HoughTransformAvg(std::unique_ptr<double[]>& arr, int arrsize, int eventnumber){


  std::cout << " Hello I am here " << std::endl;
  float m_d0_range = 120;
  float m_qOverPt_range = 0.002;
  int m_imageSize_x = 216; // i.e. number of bins in d0 (nominal 216)
  int m_imageSize_y = 216; // i.e. number of bins in q/pT (nominal 216)
  double m_step_x = (2*m_d0_range) / m_imageSize_x; // helpers (accumulator granularity)
  double m_step_y = (2*m_qOverPt_range) / m_imageSize_y;
  //std::cout << "m_step_x is " << m_step_x << " and m_step_y is " << m_step_y << endl;
  bool m_continuous = true; // assure that there is continuity of the line (i.e. middle bins in d0 are filled when one q/pT step would result in a hole)

  // if (event>2) continue;
  vector2D<std::pair<int, hit>> image(m_imageSize_x, m_imageSize_y);

  //int sizeout = sizeof(arr);
  int sizeout = arrsize/14;
  // sizeout /= 9;
  std::cout << " sizeout: " << sizeout<<  std::endl;

  std::string plotname = "outputfiles/";
  plotname.append(std::to_string(eventnumber));
  plotname.append("_outputvec_nopileup_withOuterLayers.txt");
  char nplot[plotname.size()+1];
  strcpy(nplot,plotname.c_str());
  ofstream myfile;
  myfile.open(nplot);
  //std::cout << nplot << " has been opened." << endl;
  std::vector<std::vector<int>> pix_list(m_imageSize_x*m_imageSize_y);
  int pixel = 0;

  for(int ihit1=0; ihit1<sizeout; ihit1++){
    //for(int ihit2=ihit1+1; ihit2<sizeout; ihit2++){
    for(int ihit2=1; ihit2<sizeout; ihit2++){

      const pvec p1 {{arr[14*ihit1+4], arr[14*ihit1+5]}}; // x and y for hit1
      std::cout << "x for hit 1 of " << arr[14*ihit1+1] << " is " << arr[14*ihit1+4] << " and y is " << arr[14*ihit1+5] << endl;
      const pvec p2 {{arr[14*ihit2+4], arr[14*ihit2+5]}};// x and y for hit1
      std::cout << "x for hit 2 of " << arr[14*ihit2+1] << " is " << arr[14*ihit2+4] << " and y is " << arr[14*ihit2+5] << endl;
      const pvec halfDiff = (p2 - p1)*0.5;
      const fp_t halfLen = length(halfDiff);

      int xbefore = -1;

      for ( int y = 1; y < m_imageSize_y+1; y++ ) {
        const fp_t qoverpt = -1.*( (y * m_step_y) + m_step_y*0.5 - m_qOverPt_range);
	std::cout << "for y " << y << " qoverpt is " << qoverpt << endl; 
        const fp_t radius = 1.0/(0.6*qoverpt);
	std::cout << "radius is " << radius << endl;
        const fp_t scale = std::copysign( std::sqrt( std::pow(radius/halfLen, 2) - 1), radius );
	std::cout << "scale is " << scale << endl;
        const pvec rprime = rotate90(halfDiff) * scale;
        const pvec center = p1 + halfDiff + rprime;
	std::cout << "p1[0] is " << p1[0] << " and p1[1] is " << p1[1] << endl;
	std::cout << "halfDiff[0] is " << halfDiff[0] << " and halfDiff[1] is " << halfDiff[1] << endl;
	std::cout << "rprime[0] is " << rprime[0] << " and rprime[1] is " << rprime[1] << endl;
	std::cout << "center[0] is " << center[0] << " and center[1] is " << center[1] << endl;
	std::cout << "length(center) is " << length(center) << endl;
        const fp_t d0 =  (std::signbit(radius) ? -1.0 : 1.0)*(length(center) - abs(radius));
	std::cout << "d0 is " << d0 << endl; 
        int x = (d0 + m_d0_range) / m_step_x;
	//        cout << " x: " << x << endl;
        if ( 1 <= x && x < m_imageSize_x) {
          if (xbefore == -1) xbefore = x;
          if ( m_continuous ) { // fill the bins along x starting from the last one filled
            const int xmin =  (xbefore < x)? xbefore: x;
	    //std::cout << "xmin is " << xmin << endl;
            const int xmax =  (xbefore < x)? x: xbefore;
            for ( int xinterpolated = xmin; xinterpolated <= xmax; ++xinterpolated) {
              pixel = (y-1)*m_imageSize_x + xinterpolated;
	      //std::cout << "y is " << y << " x interpolated is " << xinterpolated << " pixel is " << pixel << " hit1 is " << arr[13*ihit1+1] << " hit2 is " << arr[13*ihit2+1] << " max pixlist size is " << m_imageSize_x*m_imageSize_y << endl;
              pix_list[pixel].push_back(arr[14*ihit1+1]);
              pix_list[pixel].push_back(arr[14*ihit2+1]);
	      //std::cout << "added to pixlist" << endl;
              image(xinterpolated, y).first++;
	      myfile << xinterpolated << "," << y-1 << "," << image(xinterpolated,y).first << "," << arr[14*ihit1+1] << "," << arr[14*ihit2+1];
	      myfile << endl;
              // image(xinterpolated, y).second =  hits[event] ;
              // image(xinterpolated, y).second =  hits[event].x[ihit2] ;
              pixel=0;
            }
          } else {
	    //std::cout << "x is " << x << endl;
	    pixel = (y-1)*m_imageSize_x + x;
	    pix_list[pixel].push_back(ihit1);
	    pix_list[pixel].push_back(ihit2);
	    image(x, y).first++;
	    myfile << x << "," << y-1 << "," << image(x,y).first << "," << arr[14*ihit1+1] << "," << arr[14*ihit2+1];
	    myfile << endl;
	    // image(x, y).second =  hits[event] ;
	    // image(x, y).second.insert( hits[event].x );
	    // image(x, y).second.insert( hits[event].x );
          }
          xbefore = x;
        }
      }
    }
  }

  //sort into four categories:
  //1. Events with 1 hit per layer (Type 0)
  //2. Events with missing layers only (Type 1)
  //3. Events with extra hits only (Type 2)
  //4. Events with both missing layers and extra hits (Type 3)

  int Type = 0;
  double layerarr[sizeout];
  for(int hitties=0; hitties<sizeout; hitties++){
    layerarr[hitties] = arr[14*hitties+2];
  }
  std::sort(layerarr,layerarr+sizeout);
  if(layerarr[0] != 0){Type = 1;}
  for(int i=1; i<sizeout; i++){
    if(layerarr[i-1]+1 == layerarr[i]){
      continue;
    }
    if((layerarr[i-1]+1 > layerarr[i]) && Type == 0){
      Type = 2;
    }
    if((layerarr[i-1]+1 > layerarr[i]) && Type == 1){
      Type = 3;
    }
    if((layerarr[i-1]+1 < layerarr[i]) && Type == 0){
      Type = 1;
    }
    if((layerarr[i-1]+1 < layerarr[i]) && Type == 2){
      Type = 3;
    }
    if(Type == 3) break;
  }
  //std::cout << "Type is " << Type << endl;  
  int nroads = 0;
  int nb50 = 0;
  ofstream roadfile;
  roadfile.open("outputfiles/outputroad_nopileup.txt",std::ios_base::app);
  for (int y = 0; y < m_imageSize_y; y++) {
    for (int x = 0; x < m_imageSize_x; x++) {
      float nmatchedhits = 0;
      int pix = y*m_imageSize_x + x;
      //if (passThreshold(image, x, y, m_step_x, m_d0_range, sizeout) && isLocalMaxima( image, x, y, m_imageSize_x, m_imageSize_y) ) {
      if (passThreshold(image, x, y, m_step_x, m_d0_range, sizeout)){
	pvec testresult = LocalMaximaTest(image,x,y,m_imageSize_x,m_imageSize_y,m_step_x,m_step_y,m_d0_range,m_qOverPt_range);
	//std::cout << "testresult[x] is " << testresult[0] << endl;
	if(testresult[0] != -1){
	  nroads++;
	  std::cout << "testresult[x] is " << testresult[0] << " and testresult[y] is " << testresult[1] << endl;
	  cout << " d0: " << testresult[0] << " truthd0: " << arr[11]  << " resolution d0 :" << (arr[11] - testresult[0] ) << " q/pt " << testresult[1] << " truth q/pT: " << arr[9] / arr[10]<< " resolution q/pT :" << (arr[9] / arr[10] - testresult[1]) <<  endl;
	  for(int i=0; i<pix_list[pix].size(); i++){
	    std::cout << pix_list[pix][i] << ",";
	  }
	  for(int j=0; j<sizeout; j++){
	    for(int i=0; i<pix_list[pix].size(); i++){
	      if(j == pix_list[pix][i]){
		if(arr[14*j+7] == arr[14*j+8]){ //if hbarcode matches track barcode
		  nmatchedhits++;
		}
	      }
	    }
	  }
	  std::cout << "Barcode fraction for road is " << nmatchedhits/pix_list[pix].size() << endl;
	  if(nmatchedhits/pix_list[pix].size() > 0.5){
	    nb50++;
	  } 
	}
      }
    }
  }

  ofstream myfile2;
  myfile2.open("outputfiles/outputarr_nopileup_withOuterLayers.txt",std::ios_base::app);

  for(int i=0; i<sizeout; i++){
    arr[14*i+13]=nroads;
    //std::cout << "Efficiency for hit " << i << " is " << arr[14*i+12]/arr[14*i+13] << endl;
    myfile2 << arr[14*i] << " " << arr[14*i+1] << " " << arr[14*i+2] << " " << arr[14*i+3] << " " << arr[14*i+4] << " " << arr[14*i+5] << " " << arr[14*i+6] << " " << arr[14*i+7] << " " << arr[14*i+8] << " " << arr[14*i+9] << " " << arr[14*i+10] << " " << arr[14*i+11] << " " << nb50 << " " << arr[14*i+13] << " " << Type << endl;       }

}

#endif
