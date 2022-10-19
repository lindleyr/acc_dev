#include <dirent.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include "HoughHelper.h"
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
void GetInfoFromFiles(string hitsFile, string partsFile, hit *hits, particle *particles){
  // initialize structs (maybe cleaner way?)
  cout << " i have started initialize " << endl;

  for (int i = 0; i < 10000; ++i){
    for (int j = 0; j < 100; ++j){
      hits[i].layer[j] = -9999; hits[i].r[j] = -9999; hits[i].x[j] = -9999; hits[i].y[j] = -9999; hits[i].z[j] = -9999;
      particles[i].barcode[j] = -9999; particles[i].charge[j] = -9999; particles[i].pt[j] = -9999; particles[i].d0[j] = -9999;
    }
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
    event1_part = event_part;
  }
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
    int m_threshold = 8;
    int m_threshold50 = 8;
    if ( std::abs(d0) < 50.0 && count >= m_threshold50 ) return true;
    if ( std::abs(d0) >= 50.0 && count >= m_threshold ) return true;

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

void SelectEvents(hit *hits, particle *particles){

  const double m_acceptedDistanceBetweenLayersMin = 200; // min R disstance for hits pair filtering
  const double m_acceptedDistanceBetweenLayersMax = 600;

  float m_d0_range = 120;
  float m_qOverPt_range = 0.002;
  int m_imageSize_x = 216; // i.e. number of bins in d0
  int m_imageSize_y = 216; // i.e. number of bins in q/pT
  double m_step_x = (2*m_d0_range) / m_imageSize_x; // helpers (accumulator granularity)
  double m_step_y = (2*m_qOverPt_range) / m_imageSize_y;
  bool m_continuous = true; // assure that there is continuity of the line (i.e. middle bins in d0 are filled when one q/pT step would result in a hole)

  for(int event=0; event<10000; event++){
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




#endif
