#pragma once

#pragma unmanaged 

//Eigen 
#include "./3rdparty/Eigen/Dense"
#include "./3rdparty/Eigen/Geometry"
#include "./3rdparty/Eigen/Sparse"

typedef Eigen::Vector2i EVec2i;
typedef Eigen::Vector2f EVec2f;
typedef Eigen::SparseMatrix<float> ESpMat;
typedef Eigen::Triplet<float>      ETrip;

#include <vector>
#include <iostream>
#include <algorithm>


inline void Trace(EVec2f& p)
{
	std::cout << p[0] << "" << p[1] << std::endl;
}




//実数解を一つ返す
//solve x^3 + a x^2 + b x + c
inline float SolveCubicFunction
(
	const float a,
	const float b,
	const float c
)
{
	const float p = b - a * a /  3.0f;
	const float q = 2 * a*a*a / 27.0f - a * b / 3.0f + c;

	const float tmp = q * q / 4.0f + p * p * p / 27.0f;
	if( tmp < 0) return 0;

	const float mTmp = - q / 2.0f + sqrt(tmp);
	const float nTmp = - q / 2.0f - sqrt(tmp);
	const float m = (mTmp < 0) ?  - pow( - mTmp, 1/3.0f) :  pow(mTmp, 1/3.0f);
	const float n = (nTmp < 0) ?  - pow( - nTmp, 1/3.0f) :  pow(nTmp, 1/3.0f);

	return - a / 3.0f + m + n;
}





inline float TriangleArea(
		const EVec2f &p0, 
		const EVec2f &p1, 
		const EVec2f &p2)
{
	EVec2f v1 = p1-p0;
	EVec2f v2 = p2-p0;
	const float cross = v1[0] * v2[1] - v1[1] * v2[0];
	return 0.5f * std::abs(cross);
}





// input : cps        2D control points,        size:N)
// output: out_cps    2D Bezier control points, size:3 * N, each 3 construct quad bezier)
// output: out_points 2D curve
//
// note
// 各制御点 cps[i] に対して、2次ベジエ制御点 Ci0, Ci1, Ci2を生成する
//
// ただし，
// C_i1 は上記の C1[i]に対応
// C_i0, C_i+1,0 については、『C_i0 = C_i+1,0 = (1-λi) Ci + λi Ci+1 』と計算できるのでλiとC1[i]のみ保持する
//
//https://qiita.com/Rijicho_nl/items/05ee4c8d77e99e29daa5
//
const int ITER_NUM = 15;	

inline void compute_kCurves
(
	const std::vector<EVec2f> &cps,

	int  num_sample, //sampling number of each curve segment
	std::vector<EVec2f> &out_cps,
	std::vector<EVec2f> &out_points
)
{
	if( cps.size() < 4 ) return;
	if( num_sample < 3) num_sample = 3;

	const int N = (int)cps.size();
	std::vector<float> lambda(N, 0.5f), ti(N);
	std::vector<EVec2f> Ci0(N), Ci1 = cps, Ci2(N);

	for( int iter = 0; iter < ITER_NUM; ++iter)
	{
		//1. update lambda (use 0.5 for the first iteration )
		// iter < ITER_NUM / 2 は、安定化のためのtrick
		if( iter != 0 && iter < ITER_NUM / 2)
		{
			for( int i = 0; i < N; ++i)
			{
				int next_i = (i+1) % N;
				float A = sqrt( TriangleArea(Ci0[i], Ci1[i     ], Ci1[next_i] ) );
				float B = sqrt( TriangleArea(Ci0[i], Ci1[i     ], Ci1[next_i] ) );
				float C = sqrt( TriangleArea(Ci1[i], Ci1[next_i], Ci2[next_i] ) );
				lambda[i] = A / (B + C);
			}
		}

		//2. update Ci2, Ci+10
		for( int i = 0; i < N; ++i)
		{
			int next_i = (i+1)%N;
			Ci2[i] = Ci0[next_i] = (1 - lambda[i]) * Ci1[  i   ]  +  
				                          lambda[i]  * Ci1[next_i];
		}

		//3. update ti
		for( int i = 0; i < N; ++i)
		{
			//solve eq(4)
			EVec2f Ci2_Ci0 = Ci2[i] - Ci0[i];
			EVec2f Ci0_pi  = Ci0[i] - cps[i];
			float a = Ci2_Ci0.squaredNorm();
			float b = 3 * Ci2_Ci0.dot( Ci0_pi );
			float c = (3 * Ci0[i] - 2 * cps[i] - Ci2[i]).dot( Ci0_pi );
			float d = - Ci0_pi.squaredNorm();

			ti[i] = (a==0) ? 0.5f : SolveCubicFunction(b/a, c/a, d/a);
		}

		//4. update C1
		ESpMat A(N,N);
		Eigen::VectorXf b1(N), b2(N);
		std::vector< Eigen::Triplet<float> > entries; //{row, col, val}

		for( int i = 0; i < N; ++i)
		{
			int nextI = (i+1) % N;
			int prevI = (i-1) < 0 ? N - 1 : i - 1;
			float t = ti[i];
			float a = (1-lambda[prevI]) * (1 - t) * (1 - t);
			float b = lambda[i] * t * t;
			float c = lambda[prevI] * (1 - t) * (1 - t) + (2-(1+lambda[i]) * t) * t;

			entries.push_back( ETrip(i, prevI, a) );
			entries.push_back( ETrip(i, i    , c) );
			entries.push_back( ETrip(i, nextI, b) );
			b1[i] = cps[i][0];
			b2[i] = cps[i][1];
		}

		A.setFromTriplets(entries.begin(), entries.end());
		// solve Ax = b
		Eigen::SparseLU<ESpMat> LU(A);  
		Eigen::VectorXf x = LU.solve(b1);
		Eigen::VectorXf y = LU.solve(b2);

		for( int i = 0; i < N; ++i) Ci1[i] << x[i], y[i];

	}

	out_cps.clear();
	for( int i = 0; i < N; ++i)
	{
		out_cps.push_back( Ci0[i] );
		out_cps.push_back( Ci1[i] );
		out_cps.push_back( Ci2[i] );
	}

	out_points.clear();
	for( int i = 0; i < N; ++i)
	{
		const EVec2f &C0 = Ci0[i];
		const EVec2f &C1 = Ci1[i];
		const EVec2f &C2 = Ci2[i];

		for( int j=0; j < num_sample; ++j)
		{
			float t = j * 1.0f / num_sample;
			EVec2f p = (1 - t) * (1 - t) * C0 + 2 * t * (1 - t) * C1 + t * t * C2;
			out_points.push_back( p );
		}
	}
}

//memo 
// 
// 参照論文: "Kappa-Curves: Interpolating B-splines with Local Curvature Control" (2017)
// 
// Step 4の線形システムについて
// 
// 論文中の式 pi = α c(i-1,1) + β c(i,1) + γ c(i+1,1) は以下のように導出される。
//
// ★ 2次元ベジエの式3より
// pi = (1 - ti)^2 * c(i,0) + 2 ti (1 - ti) * c(i,1) + ti^2 * c(i,2)
// 
// ★接続条件の定義 (Eq.5)
// c(i,2) =     (1 - λi) * c( i ,1) +     λi * c(i+1,1)
// c(i,0) = (1 - λ(i-1)) * c(i-1,1) + λ(i-1) * c( i ,1),  c(i,0) = c(i-1,2)
//
// 上記 Eq.5 を Eq.3 へ代入して整理 
// pi = [ (1-λ(i-1))(1-ti)^2 ]                       * c(i-1,1)  <-- α(i)
//    + [ λ(i-1)(1-ti)^2 + 2ti(1-ti) + (1-λi)ti^2 ] * c( i ,1)  <-- β(i)
//    + [ λi * ti^2 ]                                * c(i+1,1)  <-- γ(i)
//
//
// ◯ 開いた曲線の場合 p0, p1, ..., pN に対して、
// (1)p0, pN は端点としてよけ、p1～p(N-1)のそれぞれに2Dベジエ曲線を作る
// (2) Step 4の線形システムは、以下のようになる
// 
// 
// (*)p1 において  c(1,0),c(1,1), c(1,2)によりベジエができるが、c(1,0)が端点（p0) になるので、、、
// 
// pi = (1 - ti)^2 * c(i,0) + 2 ti (1 - ti) * c(i,1) + ti^2 * c(i,2)
// c(i,2) =     (1 - λi) * c( i ,1) +     λi * c(i+1,1)
// c(i,0) = p0
// これを代入すると、、 
// pi = (1 - ti)^2 * p0 
//    + [ 2ti(1 - ti) + (1 - λi)ti^2 ] * c( i ,1)  
//    + [ λi * ti^2 ]                  * c(i+1,1)
// 
// (*)p(N-1)では c(N-1,0),c(N-1,1), c(N-1,2)によりベジエができるが、c(N-1,2)が端点（pN) になるので、、、
// 
// pi = (1 - ti)^2 * c(i,0) + 2 ti (1 - ti) * c(i,1) + ti^2 * c(i,2)
// c(i,0) = (1 - λ(i-1)) * c(i-1,1) + λ(i-1) * c( i ,1)
// c(i,2) = pN
// これを代入すると、、 
// pi = [ (1 - λ(i-1))(1 - ti)^2 ]         * c(i-1,1)
//    + [ λ(i-1)(1 - ti)^2 + 2ti(1 - ti) ] * c( i ,1)
//    + [ ti^2 ]                            * pN
// 
// N=1の時の対応
// pi = (1 - ti)^2 * start_cp + 2 ti (1 - ti) * c(1,1) + ti^2 * end_cp
// c(1,1) = ( pi - (1 - ti)^2 * start_cp - ti^2 * end_cp ) / ( 2 ti (1 - ti) )



//
// input : cps        2D control points,        size: N 
// output: out_cps    2D Bezier control points, size: 3 * (N-2), each 3 construct quad bezier
// output: out_points 2D curve
//
// note
// 各制御点 cps[i] に対して、2次ベジエ制御点 Ci0, Ci1, Ci2を生成する
//
// ただし，
// cps[1  ]に対する曲線は、cps[0]=C10, C11, C12で決定される
// cps[N-2]に対する曲線は、CN-2,0, CN-2,1, CN-2,2=cps[N-1]で決定される


//TODO N=3のと起動するか？
//TODO 置くまで届いていない理由が不明

inline void compute_kCurves_open
(
	const std::vector<EVec2f>& _cps,

	int  num_sample, //sampling number of each curve segment
	std::vector<EVec2f>& out_cps,
	std::vector<EVec2f>& out_points
)
{
	if (_cps.size() < 3) return;
	if (num_sample < 3) num_sample = 3;

	//N 本の 2次ベジエを生成する
	const int N = (int)_cps.size() - 2;
  const EVec2f start_cp = _cps.front();
  const EVec2f end_cp   = _cps.back();

  std::vector<EVec2f> cps(_cps.begin() + 1, _cps.end() - 1);
	std::vector<float>  lambda(N, 0.5f), ti(N);
	std::vector<EVec2f> Ci0(N), Ci1 = cps, Ci2(N);
	
  Ci0[0] = start_cp;
  Ci2[N - 1] = end_cp;

	for (int iter = 0; iter < ITER_NUM; ++iter)
	{
		//1. update lambda (use 0.5 for the first iteration )
		// iter < ITER_NUM / 2 は、安定化のためのtrick
		// lambda[N-1]は不要なことに注意
		if (iter != 0 && iter < ITER_NUM / 2)
		{
			for (int i = 0; i < N-1; ++i)
			{
				float A = sqrt(TriangleArea(Ci0[i], Ci1[i]    , Ci1[i + 1]));
				float B = sqrt(TriangleArea(Ci0[i], Ci1[i]    , Ci1[i + 1]));
				float C = sqrt(TriangleArea(Ci1[i], Ci1[i + 1], Ci2[i + 1]));
				lambda[i] = A / (B + C);
			}
		}

		//2. update Ci2, Ci+1,0
		for (int i = 0; i < N-1; ++i)
		{
			Ci2[i] = Ci0[i + 1] = (1 - lambda[i]) * Ci1[i] + lambda[i] * Ci1[i + 1];
		}
		Ci0[  0  ] = start_cp;
		Ci2[N - 1] = end_cp;

		//3. update ti
		for (int i = 0; i < N; ++i)
		{
			//solve eq(4)
			EVec2f Ci2_Ci0 = Ci2[i] - Ci0[i];
			EVec2f Ci0_pi  = Ci0[i] - cps[i];
			float a = Ci2_Ci0.squaredNorm();
			float b = 3 * Ci2_Ci0.dot(Ci0_pi);
			float c = (3 * Ci0[i] - 2 * cps[i] - Ci2[i]).dot(Ci0_pi);
			float d = -Ci0_pi.squaredNorm();
			ti[i] = (a == 0) ? 0.5f : SolveCubicFunction(b / a, c / a, d / a);
			if (std::isnan(ti[i])) ti[i] = 0.5f;
			ti[i] = std::max(1e-4f, std::min(1.0f - 1e-4f, ti[i]));
		}

		//4. update C1 
		ESpMat A(N, N);
		Eigen::VectorXf b1(N), b2(N);
		std::vector< Eigen::Triplet<float> > entries; //{row, col, val}

		for (int i = 0; i < N; ++i)
		{
			const float t = ti[i];
			if (i ==0 && i == N-1) // N =1 のとき
			{
				entries.push_back(ETrip(i, i, 2 * t * (1 - t)));
				b1[i] = cps[i][0] - (1 - t) * (1 - t) * start_cp[0] - t * t * end_cp[0];
				b2[i] = cps[i][1] - (1 - t) * (1 - t) * start_cp[1] - t * t * end_cp[1];
      }
			else if (i == 0)
			{
				const float b = 2 * t * (1 - t) + (1 - lambda[i]) * t * t;
				const float c = lambda[i] * t * t;
				entries.push_back(ETrip(i, i    , b));
				entries.push_back(ETrip(i, i + 1, c));
				b1[i] = cps[i][0] - (1 - t) * (1 - t) * start_cp[0];
				b2[i] = cps[i][1] - (1 - t) * (1 - t) * start_cp[1];
      }
      else if (i == N - 1)
			{
				const float a = (1 - lambda[i - 1]) * (1 - t) * (1 - t) ;
				const float b = lambda[i - 1] * (1 - t) * (1 - t) + 2 * t * (1 - t);
				entries.push_back(ETrip(i, i-1, a));
				entries.push_back(ETrip(i, i  , b));
				b1[i] = cps[i][0] - t * t * end_cp[0];
				b2[i] = cps[i][1] - t * t * end_cp[1];
			}
			else
			{
				const float a = (1 - lambda[i - 1]) * (1 - t) * (1 - t);
				const float b = lambda[i - 1] * (1 - t) * (1 - t) + (2 - (1 + lambda[i]) * t) * t;
				const float c = lambda[i] * t * t;
				entries.push_back(ETrip(i, i-1, a));
				entries.push_back(ETrip(i,  i , b));
				entries.push_back(ETrip(i, i+1, c));
				b1[i] = cps[i][0];
				b2[i] = cps[i][1];
			}
		}

		A.setFromTriplets(entries.begin(), entries.end());
		// solve Ax = b
		Eigen::SparseLU<ESpMat> LU(A);
		Eigen::VectorXf x = LU.solve(b1);
		Eigen::VectorXf y = LU.solve(b2);

		for (int i = 0; i < N; ++i) Ci1[i] << x[i], y[i];
	}

	out_cps.clear();
	for (int i = 0; i < N; ++i)
	{
		out_cps.push_back(Ci0[i]);
		out_cps.push_back(Ci1[i]);
		out_cps.push_back(Ci2[i]);
	}

	out_points.clear();
	for (int i = 0; i < N; ++i)
	{
		const EVec2f& C0 = Ci0[i];
		const EVec2f& C1 = Ci1[i];
		const EVec2f& C2 = Ci2[i];

		for (int j = 0; j < num_sample; ++j)
		{
			float t = j * 1.0f / (num_sample-1);
			EVec2f p = (1 - t) * (1 - t) * C0 + 2 * t * (1 - t) * C1 + t * t * C2;
			out_points.push_back(p);
		}
	}
}


#pragma managed 



