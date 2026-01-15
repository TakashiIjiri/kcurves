#pragma once

#pragma unmanaged 

#include "./3rdparty/Eigen/Dense"
#include "./3rdparty/Eigen/Geometry"

typedef Eigen::Vector2i EVec2i;
typedef Eigen::Vector2f EVec2f;

#include <vector>
#include <iostream>
#include <algorithm>

namespace catmullromspline
{

  //https://ja.wikipedia.org/wiki/Catmull-Rom%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3%E6%9B%B2%E7%B7%9A
  //https://ja.wikipedia.org/wiki/Centripetal_Catmull-Rom%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3
  //https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline?utm_source=chatgpt.com

  // alpha = 0.0 -> uniform
  // alpha = 0.5 -> centripetal
  std::vector<EVec2f> compute_catmullrom_spline
  (
    const std::vector<EVec2f>& cps,
    float alpha = 0.5,
    int num_samples_per_segment = 30
  )
  {
    const int N = static_cast<int>(cps.size());
    
    num_samples_per_segment = std::max(4, num_samples_per_segment);

    std::vector<EVec2f> points;

    if ( N < 2 ) return points;
    
    if ( N == 2)
    {
      for (int j = 0; j < num_samples_per_segment; ++j)
      {
        float t = (float) j / (num_samples_per_segment - 1);
        EVec2f point = (1.0f - t) * cps[0] + t * cps[1];
        points.push_back(point);
      }
      return points;
    }

    const float EPS = 1e-5f;
    const EVec2f start_cp = cps[0  ] + (cps[ 0 ] - cps[ 1 ]);
    const EVec2f end_cp   = cps[N-1] + (cps[N-1] - cps[N-2]);

    for (int i = 0; i < N - 1; ++i)
    {
      const EVec2f &p0 = (i ==  0 ) ? start_cp : cps[i - 1];
      const EVec2f &p1 = cps[i    ];
      const EVec2f &p2 = cps[i + 1];
      const EVec2f &p3 = (i == N-2) ? end_cp   : cps[i + 2];
      const float t0 = 0.0f;
      const float t1 = t0 + std::pow(std::max((p1 - p0).norm(), EPS), alpha);
      const float t2 = t1 + std::pow(std::max((p2 - p1).norm(), EPS), alpha);
      const float t3 = t2 + std::pow(std::max((p3 - p2).norm(), EPS), alpha);

      for (int j = 0; j < num_samples_per_segment; ++j)
      {
        float t = t1 + ((float)j / (float)num_samples_per_segment) * (t2 - t1);
        EVec2f A1 = ((t1 - t) / (t1 - t0)) * p0 + ((t - t0) / (t1 - t0)) * p1;
        EVec2f A2 = ((t2 - t) / (t2 - t1)) * p1 + ((t - t1) / (t2 - t1)) * p2;
        EVec2f A3 = ((t3 - t) / (t3 - t2)) * p2 + ((t - t2) / (t3 - t2)) * p3;

        EVec2f B1 = ((t2 - t) / (t2 - t0)) * A1 + ((t - t0) / (t2 - t0)) * A2;
        EVec2f B2 = ((t3 - t) / (t3 - t1)) * A2 + ((t - t1) / (t3 - t1)) * A3;  

        EVec2f p = ((t2 - t) / (t2 - t1)) * B1 + ((t - t1) / (t2 - t1)) * B2;
        points.push_back(p);  
      }
    }
    points.push_back(cps.back());
    return points;

  }








};


#pragma managed 



